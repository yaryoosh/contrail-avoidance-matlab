%  CONTRAIL AVOIDANCE OPTIMISER  —  London Heathrow (LHR) -> New York (JFK)

%  1. Loads ERA5 pressure-level data (t, q, u, v) from a NetCDF file
%  2. Builds a 3D map of persistent-contrail regions using the
%     Schmidt–Appleman criterion + ice supersaturation (RHi >= 100%)
%  3. Builds a 3D graph of possible routes (along-track x lateral x level)
%  4. Costs every segment with a physics-based fuel model (drag polar,
%     TSFC, real winds, climb penalty) and contrail distance
%  5. Runs Dijkstra for a sweep of trade-off weights (lambda) to get a
%     fuel vs contrail Pareto front
%  6. Produces: contrail map, 3D route plot, Pareto plot, altitude
%     profile, and a rotating GIF for LinkedIn

clear; clc; close all;

%% User settings
dataFile  = 'contrails_data_sheet.nc';  
timeIndex = 1;              

% Route end points [lat, lon] in degrees
orig = [51.4700,  -0.4543];  origName = 'LHR';
dest = [40.6413, -73.7781];  destName = 'JFK';

% Route graph resolution
nAlong        = 61;                    % points along the great circle
nOffSide      = 6;                     % lateral offsets each side of track
offStep_km    = 40;                    % lateral spacing between offsets (km)
cruiseLev_hPa = [200 225 250 300];     % allowed cruise levels (~FL300-FL390)
refLev_hPa    = 250;                   % level for the "fixed-level" reference

% Aircraft model (A350-900-like, mid-cruise mass)
ac.m        = 230000;           % mass (kg)
ac.S        = 442;              % wing area (m^2)
ac.CD0      = 0.018;            % zero-lift drag coefficient
ac.k        = 1/(pi*9.49*0.85); % induced drag factor 1/(pi*AR*e)
ac.Mach     = 0.85;             % cruise Mach number
ac.TSFC     = 1.47e-5;          % thrust-specific fuel consumption (kg/(N s))
ac.etaClimb = 0.30;             % fuel energy -> potential energy efficiency

% Contrail physics
eta_SAC  = 0.35;   % overall propulsion efficiency used in Schmidt–Appleman
RHiScale = 1.00;   % >1 (e.g. 1.10) to correct ERA5's known dry bias

% Trade-off sweep: lambda = kg of fuel you'd "pay" to avoid 1 km of contrail
lambdas = [0, logspace(-1, 3, 30)];

% Output
saveFigures = true;
gifName     = 'contrail_routes.gif';

%% Constants
c.g = 9.80665;  c.R = 287.05;  c.gam = 1.4;  c.Q = 43.2e6;

%% =Loading the ERA5 file / Contrails_data_sheet file
fprintf('Loading %s ...\n', dataFile);
lon  = readCoord(dataFile, {'longitude','lon'});
lat  = readCoord(dataFile, {'latitude','lat'});
plev = readCoord(dataFile, {'pressure_level','level','isobaricInhPa'});
if max(plev) > 2000, plev = plev/100; end          % Pa -> hPa if needed

T = readField(dataFile, {'t','temperature'},           timeIndex); % K
q = readField(dataFile, {'q','specific_humidity'},     timeIndex); % kg/kg
u = readField(dataFile, {'u','u_component_of_wind'},   timeIndex); % m/s
v = readField(dataFile, {'v','v_component_of_wind'},   timeIndex); % m/s

% Make all grid vectors ascending (ERA5 stores latitude north -> south)
lon(lon > 180) = lon(lon > 180) - 360;
[lon,  iLo] = sort(lon);
[lat,  iLa] = sort(lat);
[plev, iLv] = sort(plev);
T = T(iLo,iLa,iLv);  q = q(iLo,iLa,iLv);
u = u(iLo,iLa,iLv);  v = v(iLo,iLa,iLv);

fprintf('Grid: %d lon x %d lat x %d levels\n', numel(lon), numel(lat), numel(plev));
fprintf('Levels (hPa): %s\n', mat2str(plev'));
fprintf('Temperature range: %.1f to %.1f K\n', min(T(:)), max(T(:)));

%% Persistent control mask
P = reshape(plev*100, 1, 1, []);                 % pressure in Pa, [1 1 nLev]
[mask, RHi, TLM] = contrailMask(T, q, P, eta_SAC, RHiScale);

fprintf('\nLevel   T_LM(degC)   %%ice-supersat.   %%persistent contrail\n');
for k = 1:numel(plev)
    fprintf('%4d hPa   %7.1f        %6.1f            %6.1f\n', plev(k), ...
        TLM(1,1,k)-273.15, 100*mean(RHi(:,:,k) >= 1, 'all'), ...
        100*mean(mask(:,:,k), 'all'));
end
if ~any(mask(:))
    warning('No persistent-contrail regions found. Try another day, or RHiScale = 1.1.');
end

% Interpolants (lon, lat, pressure in Pa)
gv  = {lon, lat, plev*100};
F.T = griddedInterpolant(gv, T,            'linear', 'nearest');
F.U = griddedInterpolant(gv, u,            'linear', 'nearest');
F.V = griddedInterpolant(gv, v,            'linear', 'nearest');
F.M = griddedInterpolant(gv, double(mask), 'linear', 'nearest');

%% Build route grid
cruiseLev_hPa = sort(cruiseLev_hPa(cruiseLev_hPa >= min(plev) & ...
                                   cruiseLev_hPa <= max(plev)));
if numel(cruiseLev_hPa) < 2, error('Need at least 2 cruise levels inside the data.'); end
pLev = cruiseLev_hPa(:)*100;          % Pa, column
nLev = numel(pLev);

[latC, lonC] = gcPoints(orig(1), orig(2), dest(1), dest(2), nAlong);
brg = bearingDeg(latC(1:end-1), lonC(1:end-1), latC(2:end), lonC(2:end));
brg = [brg; brg(end)];

offs = (-nOffSide:nOffSide) * offStep_km * 1e3;   % metres, row
nOff = numel(offs);
j0   = nOffSide + 1;                               % centre-line index
[nodeLat, nodeLon] = destPoint(latC, lonC, brg + 90, offs);   % nAlong x nOff

nN  = nAlong*nOff*nLev;
id  = @(i,j,l) sub2ind([nAlong nOff nLev], i, j, l);
src = nN + 1;  snk = nN + 2;

% Edges: each step forward may move 1 offset sideways and/or 1 level up/down
[I,J,L,DJ,DL] = ndgrid(1:nAlong-1, 1:nOff, 1:nLev, -1:1, -1:1);
J2 = J + DJ;  L2 = L + DL;
ok = J2>=1 & J2<=nOff & L2>=1 & L2<=nLev;
I = I(ok); J = J(ok); L = L(ok); J2 = J2(ok); L2 = L2(ok);

s = id(I, J, L);
t = id(I+1, J2, L2);
latA = nodeLat(sub2ind([nAlong nOff], I,   J ));
lonA = nodeLon(sub2ind([nAlong nOff], I,   J ));
latB = nodeLat(sub2ind([nAlong nOff], I+1, J2));
lonB = nodeLon(sub2ind([nAlong nOff], I+1, J2));

fprintf('\nCosting %d route segments ...\n', numel(s));
[fuelE, conE, distE] = segCost(latA, lonA, pLev(L), latB, lonB, pLev(L2), F, ac, c);

% Virtual source/sink: route must start/end on the centre line, any level
lv  = (1:nLev)';
sS  = repmat(src, nLev, 1);            tS = id(ones(nLev,1),    j0*ones(nLev,1), lv);
sK  = id(nAlong*ones(nLev,1), j0*ones(nLev,1), lv);   tK = repmat(snk, nLev, 1);
s     = [s; sS; sK];         t     = [t; tS; tK];
fuelE = [fuelE; zeros(2*nLev,1)];
conE  = [conE;  zeros(2*nLev,1)];
distE = [distE; zeros(2*nLev,1)];

% digraph re-orders edges, so store original edge index as the weight first
nE   = numel(s);
G0   = digraph(s, t, (1:nE)', nN + 2);
perm = G0.Edges.Weight;

%% Optimise lambda sweep
res = struct('lambda',{},'fuel',{},'con',{},'dist',{},'lat',{},'lon',{},'p',{},'alt',{});
for k = 1:numel(lambdas)
    w = fuelE + lambdas(k)*conE + 1e-6;      % tiny epsilon keeps weights > 0
    G = G0;  G.Edges.Weight = w(perm);
    [Pth, ~, ep] = shortestpath(G, src, snk);
    e = perm(ep);
    nodes = Pth(2:end-1);
    [ii, jj, ll] = ind2sub([nAlong nOff nLev], nodes(:));
    res(k).lambda = lambdas(k);
    res(k).fuel   = sum(fuelE(e));
    res(k).con    = sum(conE(e));
    res(k).dist   = sum(distE(e))/1e3;
    res(k).lat    = nodeLat(sub2ind([nAlong nOff], ii, jj));
    res(k).lon    = nodeLon(sub2ind([nAlong nOff], ii, jj));
    res(k).p      = pLev(ll);
    res(k).alt    = isaAlt(pLev(ll))/1e3;
end

% Reference: plain great circle at a single fixed level
[~, lRef] = min(abs(cruiseLev_hPa - refLev_hPa));
pR = pLev(lRef)*ones(nAlong-1,1);
[fR, cR] = segCost(latC(1:end-1), lonC(1:end-1), pR, latC(2:end), lonC(2:end), pR, F, ac, c);
ref.fuel = sum(fR);  ref.con = sum(cR);

fuel = [res.fuel];  con = [res.con];
fuel0 = fuel(1);    con0 = con(1);
dF = 100*(fuel - fuel0)/fuel0;
if con0 > 0
    dC = 100*(con0 - con)/con0;
else
    dC = zeros(size(con));
    warning('Fuel-optimal route already has no contrails on this day — pick another day.');
end

% "Best" trade-off = cheapest route achieving >= 90% of the max reduction
kBest = find(dC >= 0.9*max(dC), 1);
if isempty(kBest) || max(dC) == 0, kBest = numel(res); end

fprintf('\n================ RESULTS ================\n');
fprintf('Fixed %d hPa great circle : fuel %7.0f kg | contrail %6.0f km\n', ...
        cruiseLev_hPa(lRef), ref.fuel, ref.con);
fprintf('Fuel-optimal (lambda = 0) : fuel %7.0f kg | contrail %6.0f km\n', fuel0, con0);
fprintf('Chosen trade-off          : fuel %7.0f kg | contrail %6.0f km (lambda = %.2f)\n', ...
        fuel(kBest), con(kBest), lambdas(kBest));
fprintf('  -> %.1f%% contrail distance avoided for %.2f%% extra fuel\n', dC(kBest), dF(kBest));
fprintf('\n lambda     fuel(kg)  contrail(km)  +fuel(%%)  -contrail(%%)\n');
for k = 1:numel(res)
    fprintf('%8.2f  %9.0f   %9.0f   %7.2f   %9.1f\n', lambdas(k), fuel(k), con(k), dF(k), dC(k));
end

save('contrail_results.mat', 'res', 'ref', 'lambdas', 'dF', 'dC', 'kBest');

%% Plots
colOpt = [0.85 0.20 0.20];   % fuel-optimal route
colAvd = [0.10 0.65 0.30];   % contrail-avoiding route
rOpt = res(1);  rAvd = res(kBest);
coast = tryCoastlines();

% Figure 1: contrail map at the level the route mostly flies
[~, iMap] = min(abs(plev - mode(rAvd.p)/100));
f1 = figure('Color','w','Position',[50 50 1100 600]); lightTheme(f1);
pcolor(lon, lat, double(mask(:,:,iMap))'); shading flat; hold on
colormap(f1, [1 1 1; 0.55 0.75 1]); clim([0 1]);
hC = patch(NaN, NaN, [0.55 0.75 1], 'EdgeColor', 'none');   % legend swatch
if ~isempty(coast), plot(coast.lon, coast.lat, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.8); end
hG = plot(lonC, latC, 'k--', 'LineWidth', 1);
hO = plot(rOpt.lon, rOpt.lat, '-', 'Color', colOpt, 'LineWidth', 5);   % thick underneath
hA = plot(rAvd.lon, rAvd.lat, '-', 'Color', colAvd, 'LineWidth', 2);   % thin on top
plot([orig(2) dest(2)], [orig(1) dest(1)], 'ko', 'MarkerFaceColor', 'k');
text(orig(2)+1, orig(1)+1, origName, 'FontWeight', 'bold');
text(dest(2)+1, dest(1)+1, destName, 'FontWeight', 'bold');
xlim([min(lon) max(lon)]); ylim([min(lat) max(lat)]);
xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
title(sprintf('Persistent-contrail regions at %d hPa and optimised routes', plev(iMap)));
legend([hC hG hO hA], {'Contrail region','Great circle','Fuel-optimal','Contrail-avoiding'}, ...
       'Location','southwest');
box on; set(gca, 'FontSize', 12, 'Layer', 'top');

% Figure 2: 3D view
altLev = isaAlt(plev*100)/1e3;
[altS, ord] = sort(altLev);
Vm = permute(double(mask(:,:,ord)), [2 1 3]);      % [lat lon alt] for meshgrid
[X, Y, Z] = meshgrid(lon, lat, altS);
f2 = figure('Color','w','Position',[80 80 1100 650]); lightTheme(f2); hold on
if any(Vm(:))
    fv = isosurface(X, Y, Z, smooth3(Vm, 'box', 3), 0.5);
    patch(fv, 'FaceColor', [0.35 0.6 1], 'EdgeColor', 'none', 'FaceAlpha', 0.35);
end
if ~isempty(coast)
    plot3(coast.lon, coast.lat, min(altS)*ones(size(coast.lon)), 'Color', [0.4 0.4 0.4]);
end
plot3(rOpt.lon, rOpt.lat, rOpt.alt, '-', 'Color', colOpt, 'LineWidth', 3);
plot3(rAvd.lon, rAvd.lat, rAvd.alt, '-', 'Color', colAvd, 'LineWidth', 3);
xlim([min(lon) max(lon)]); ylim([min(lat) max(lat)]); zlim([min(altS) max(altS)]);
xlabel('Longitude'); ylabel('Latitude'); zlabel('Altitude (km)');
title({sprintf('%s \\rightarrow %s', origName, destName), 'fuel-optimal (red) vs contrail-avoiding (green)'});
pbaspect([2.2 1.2 0.7]); view(-35, 28); grid on; box on;
camlight headlight; lighting gouraud; set(gca, 'FontSize', 12);

% Figure 3: Pareto front
f3 = figure('Color','w','Position',[120 120 750 520]); lightTheme(f3);
plot(dF, dC, '-o', 'Color', [0.2 0.2 0.2], 'MarkerFaceColor', [0.6 0.6 0.6], 'LineWidth', 1.5); hold on
plot(dF(kBest), dC(kBest), 'o', 'MarkerSize', 12, 'MarkerFaceColor', colAvd, 'MarkerEdgeColor', 'k');
text(dF(kBest), dC(kBest), sprintf('   %.0f%% fewer contrails\n   for +%.2f%% fuel', dC(kBest), dF(kBest)), ...
     'FontSize', 12, 'VerticalAlignment', 'top');
xlabel('Extra fuel burn vs fuel-optimal route (%)');
ylabel('Contrail distance avoided (%)');
title('Fuel vs contrail trade-off (Pareto front)');
grid on; box on; set(gca, 'FontSize', 12);

% Figure 4: contrail "curtain" along the track + altitude profiles
pc   = linspace(min(plev), max(plev), 80)*100;          % Pa
dAvd = alongDist(rAvd);
[PP, II] = meshgrid(pc, 1:numel(rAvd.lat));
Mc   = F.M(rAvd.lon(II), rAvd.lat(II), PP);              % contrail field under the route
f4 = figure('Color','w','Position',[150 150 1000 450]); lightTheme(f4); hold on
pcolor(dAvd, isaAlt(pc)/1e3, Mc'); shading interp
colormap(f4, [linspace(1,0.55,64)', linspace(1,0.75,64)', ones(64,1)]); caxis([0 1]);
hC4 = patch(NaN, NaN, [0.55 0.75 1], 'EdgeColor', 'none');
hA4 = plot(dAvd, rAvd.alt, '-', 'Color', colAvd, 'LineWidth', 3);
xlabel('Distance flown (km)'); ylabel('Altitude (km)');
xlim([0 max(dAvd)]); ylim([min(isaAlt(pc)) max(isaAlt(pc))]/1e3);
title('Vertical slice of contrail-forming air along the contrail-avoiding route');
legend([hC4 hA4], {'Contrail region','Contrail-avoiding route'}, 'Location', 'best');
box on; set(gca, 'FontSize', 12, 'Layer', 'top');

if saveFigures
    saveFig(f1, 'fig1_contrail_map.png');
    saveFig(f2, 'fig2_routes_3d.png');
    saveFig(f3, 'fig3_pareto.png');
    saveFig(f4, 'fig4_altitude_profile.png');
end

% Rotating GIF of the 3D figure
if isvalid(f2)
    ax2 = findobj(f2, 'Type', 'axes');
    if exist(gifName, 'file'), delete(gifName); end
    azs = linspace(-35, 325, 72);
    for k = 1:numel(azs)
        view(ax2, azs(k), 28); drawnow;
        exportgraphics(f2, gifName, 'Append', true);   % needs R2022a+
    end
    view(ax2, -35, 28);
else
    warning('3D figure was closed before the GIF step — re-run and leave it open.');
end
fprintf('\nSaved figures, %s and contrail_results.mat\n', gifName);


%Functions

function [mask, RHi, TLM] = contrailMask(T, q, P, eta, scale)
% Schmidt–Appleman criterion (Schumann 1996) + ice supersaturation.
    eps_ = 0.622; cp = 1004; Q = 43.2e6; EI = 1.25;
    G   = EI*cp.*P ./ (eps_*Q*(1 - eta));                 % Pa/K
    x   = log(G - 0.053);
    TLM = -46.46 + 9.43*x + 0.72*x.^2 + 273.15;           % K
    e   = q.*P ./ (eps_ + (1 - eps_).*q);                 % vapour pressure (Pa)
    esi = exp(9.550426 - 5723.265./T + 3.53068*log(T) - 0.00728332*T); % Murphy & Koop 2005
    RHi = scale * e ./ esi;
    mask = (T < TLM) & (RHi >= 1);
end

function [fuel, conKm, dist] = segCost(latA, lonA, pA, latB, lonB, pB, F, ac, c)
% Fuel (kg), contrail distance (km) and length (m) of each segment.
    nS = 6;  f = ((1:nS) - 0.5)/nS;                       % sample points
    latS = latA + (latB - latA).*f;
    lonS = lonA + (lonB - lonA).*f;
    pS   = pA   + (pB   - pA  ).*f;
    dist = haversine(latA, lonA, latB, lonB);
    th   = bearingDeg(latA, lonA, latB, lonB);

    Ts = F.T(lonS, latS, pS);  us = F.U(lonS, latS, pS);
    vs = F.V(lonS, latS, pS);  ms = F.M(lonS, latS, pS);

    rho = pS ./ (c.R*Ts);
    V   = ac.Mach * sqrt(c.gam*c.R*Ts);                    % true airspeed
    qd  = 0.5*rho.*V.^2;
    CL  = ac.m*c.g ./ (qd*ac.S);
    D   = qd*ac.S .* (ac.CD0 + ac.k*CL.^2);                % drag = thrust
    GS  = max(V + us.*sind(th) + vs.*cosd(th), 50);        % ground speed
    fuel = sum(ac.TSFC * D .* (dist/nS) ./ GS, 2);

    dh   = max(isaAlt(pB) - isaAlt(pA), 0);                % climb only
    fuel = fuel + ac.m*c.g*dh ./ (ac.etaClimb*c.Q);
    conKm = mean(ms, 2) .* dist/1e3;
end

function h = isaAlt(p)
% ISA pressure altitude (m) from pressure (Pa), troposphere + lower stratosphere.
    h  = 44330.77*(1 - (p/101325).^0.1902632);
    hi = p < 22632.06;
    h(hi) = 11000 + 6341.62*log(22632.06 ./ p(hi));
end

function d = haversine(lat1, lon1, lat2, lon2)
    R = 6371e3;
    a = sind((lat2-lat1)/2).^2 + cosd(lat1).*cosd(lat2).*sind((lon2-lon1)/2).^2;
    d = 2*R*asin(min(1, sqrt(a)));
end

function b = bearingDeg(lat1, lon1, lat2, lon2)
    dl = lon2 - lon1;
    y  = sind(dl).*cosd(lat2);
    x  = cosd(lat1).*sind(lat2) - sind(lat1).*cosd(lat2).*cosd(dl);
    b  = mod(atan2d(y, x), 360);
end

function [lat2, lon2] = destPoint(lat1, lon1, brg, d)
% Point reached from (lat1,lon1) going distance d (m) on bearing brg (deg).
    dd   = d/6371e3;
    lat2 = asind(sind(lat1).*cos(dd) + cosd(lat1).*sin(dd).*cosd(brg));
    lon2 = lon1 + atan2d(sind(brg).*sin(dd).*cosd(lat1), cos(dd) - sind(lat1).*sind(lat2));
end

function [lat, lon] = gcPoints(lat1, lon1, lat2, lon2, n)
% n evenly spaced points along the great circle (spherical interpolation).
    p1 = [cosd(lat1)*cosd(lon1), cosd(lat1)*sind(lon1), sind(lat1)];
    p2 = [cosd(lat2)*cosd(lon2), cosd(lat2)*sind(lon2), sind(lat2)];
    om = acos(max(-1, min(1, dot(p1, p2))));
    f  = linspace(0, 1, n)';
    Pp = (sin((1-f)*om).*p1 + sin(f*om).*p2) / sin(om);
    lat = asind(Pp(:,3));
    lon = atan2d(Pp(:,2), Pp(:,1));
end

function d = alongDist(r)
    d = [0; cumsum(haversine(r.lat(1:end-1), r.lon(1:end-1), r.lat(2:end), r.lon(2:end)))]/1e3;
end

function vn = findVar(fn, names)
    info = ncinfo(fn);
    allNames = {info.Variables.Name};
    idx = find(ismember(lower(allNames), lower(names)), 1);
    if isempty(idx)
        error('None of {%s} found in %s. Run ncdisp(''%s'') to check names.', ...
              strjoin(names, ', '), fn, fn);
    end
    vn = allNames{idx};
end

function x = readCoord(fn, names)
    x = double(ncread(fn, findVar(fn, names)));
    x = x(:);
end

function A = readField(fn, names, tIdx)
% Reads a 4D ERA5 field and returns it as [lon x lat x level] for one time.
    vn = findVar(fn, names);
    A  = double(ncread(fn, vn));
    vi = ncinfo(fn, vn);
    dn = lower({vi.Dimensions.Name});
    iLon = find(ismember(dn, {'longitude','lon'}));
    iLat = find(ismember(dn, {'latitude','lat'}));
    iLev = find(ismember(dn, {'pressure_level','level','isobaricinhpa'}));
    rest = setdiff(1:numel(dn), [iLon iLat iLev]);
    A  = permute(A, [iLon iLat iLev rest]);
    sz = size(A, 1:3);
    A  = reshape(A, sz(1), sz(2), sz(3), []);
    A  = A(:, :, :, min(tIdx, size(A, 4)));
end

function coast = tryCoastlines()
% Coastlines from Mapping Toolbox if available, otherwise skip.
    coast = [];
    try
        S = load('coastlines');
        coast.lat = S.coastlat;  coast.lon = S.coastlon;
    catch
    end
end

function lightTheme(f)
% Forces light theme (MATLAB R2025a+ dark mode makes text unreadable on white).
    try theme(f, 'light'); catch, end
end

function saveFig(f, name)
    if ~isgraphics(f)
        warning('%s skipped: the figure window was closed before saving.', name);
        return
    end
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    if exist('exportgraphics', 'file')
        exportgraphics(f, name, 'Resolution', 200);
    else
        saveas(f, name);
    end
end
