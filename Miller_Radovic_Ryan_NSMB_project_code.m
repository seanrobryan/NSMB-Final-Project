% Name: Jacob Miller, Stefan Radovic, & Sean Ryan
% Date: April 30, 2022
% Filename: Miller_Radovic_Ryan_NSMB_project_code.m
% Description: All code used to the completion of the final project in 
% BME: 5441 : Numerical and Statistical Methods for Bioengineers

clear
close all
%% Assumptions
wt_Ed = 200;  % assumed weight of Ed Sander (lbs)
wt_bike = 29;  % assumed weight of bike (lbs)
wt = wt_Ed + wt_bike;  % total weight (lbs)
m = wt*.45359237;  % kg
Gas_yield = 1.3*10^5;  % kJ/gal (citation needed)
g = 9.81;  % m/s^2

Rearth = 6371000; % radius of the Earth in meters, approximately correct for Iowa
milesPerMeter = 0.000621371; % Conversion factor meters to miles
ftsq_to_msq = 0.0929; % m^2 / ft^2

humanEnergyEff = 0.25; % Human Energy Efficiency Factor

%% Data Intake
Crr = [0.004; 0.004; 0.004; 0.008]; % Rolling resistance coeff 
                                                                            % for bike on asphault x 3, gravel x1

% drag chart values (source:
% https://www.cse.iitk.ac.in/users/amit/books/whitt-1982-bicycling-science.html)
% Col 1 := Area in ft^2, Col 2 := Drag coefficient
chartVals = [
    3.9  0.88;       % Full Racing
    4.3 1.0;        % Touring
    5.5 1.1        % European Commuter
    ];
chartVals(:, 1) = chartVals(:, 1) * ftsq_to_msq;

 % Increment between first and last known area values
areaRange = chartVals(1, 1):0.01:chartVals(end, 1);
len = length(areaRange);
linearDragVals = zeros(1, len);
lagrangeDragVals = zeros(1, len);
linearDragVals(1) = chartVals(1, 2); 
linearDragVals(end) = chartVals(end, 2);
lagrangeDragVals(1) = linearDragVals(1);

%% Interpolation
% Linear interpolation of drag (x) from known (x0, y0) and (y0, y1)
% incrementing x by 0.01
for i=2:len - 1
    linearDragVals(i) = linearDragVals(1) + (areaRange(i) - areaRange(1)) * ...
        (linearDragVals(end) - linearDragVals(1))/(areaRange(end) - areaRange(1));
end

% Lagrange Interpolation
x0 = areaRange(1);
x1 = areaRange(2);
x2 = areaRange(3);
y0 = linearDragVals(1);
y1 = chartVals(2, 2);  % y1 = linearDragVals(2);
y2 = linearDragVals(3);

for i=2:len-1
    x = areaRange(i);
    lagrangeDragVals(i) = ( ((x - x1) * (x - x2))/((x0 - x1)*(x0 - x2)) ) * y0 + ...
        ( ((x - x0) * (x - x2))/((x1 - x0)*(x1 - x2)) ) * y1 + ...
            ( ((x - x0) * (x - x1))/((x2 - x0)*(x2 - x1)) ) * y2 ;
end

rides = [ 
    "Century_Ride_with_Rolf.gpx"; 
    "Florida_Man_.gpx"; 
    "Zwift_ZHQ_FutureWorks_Crit_City_Race_Anti_sandbagging_B_.gpx";
    "Jinglecross_gravel_race.gpx"];

stravaEnergyOutputEstimates = [
    2865; 
    2516; 
    377;
    2937];  % kj

nRides = length(rides);

% initialize variables (preallocate)
linearAreaEst = zeros(1, nRides);
linearCDragEst = zeros(1, nRides);
linearEnergyExpenditureEst = zeros(1, nRides);
linearMPGEst = zeros(1, nRides);

lagrangeAreaEst = zeros(1, nRides);
lagrangeCDragEst = zeros(1, nRides);
lagrangeEnergyExpenditureEst = zeros(1, nRides);
lagrangeMPGEst = zeros(1, nRides);

linearPercentErr = zeros(1, nRides);
lagrangePercentErr = zeros(1, nRides);

%% Calculations loop
for idx=1:length(rides)
    ride = rides(idx);
    P = gpxread("data/" + ride);

    stravaEst = stravaEnergyOutputEstimates(idx);

    latitude = P.Latitude;
    longitude = P.Longitude;
    elevation = P.Elevation;
    
%     wm = webmap('World Imagery');
%     s = geoshape(latitude,longitude);
%     wmline(s,'Color','red','width',3)
    
    timeStr = strrep(P.Time, 'T', ' ');
    timeStr = strrep(timeStr, 'Z', '');
    time = datetime(timeStr);
    dtime_sec = seconds(diff(time));
    
    % Cumulative Elapsed time
    cum_Elapsed_time_sec = cumsum(dtime_sec);
    elaspse_t=seconds(cum_Elapsed_time_sec(end));
    elaspse_t.Format='hh:mm:ss';

    % Calculate distance from latitude and longitude data. 
    % Haversine formula - short for half reverse sine, where 1-cos(theta) is
    % reverse sine, and (1-cos(theta))/2  = sin^2(theta/2)
     
    phi = deg2rad(latitude);  % Could do it phi = latitude*pi/180;
    lam = deg2rad(longitude);
    dphi = diff(phi); % Could do it dphi = phi(2:end)-phi(1:end-1);
    dlam = diff(lam);   
     
    a = sin(dphi/2).^2 + cos(phi(1:end-1)).*cos(phi(2:end)).*sin(dlam/2).^2;
    ddist_xy = 2.*Rearth.*asin(sqrt(a));  
    
    delev = diff(elevation);
    ddist_xyz = sqrt(ddist_xy.^2 + delev.^2);
    
    dist_cum=cumsum(ddist_xyz);  % Cumulative total distance
    dist_traveled = dist_cum(end);  % meters
    
    %% Angle Calc.
    theta = zeros(1, length(delev));
    for i = 1:length(delev)
        if delev(i) ~= 0
            theta(i) = atan(delev(i)/ddist_xy(i));
        else
            theta(i) = atan(0);
        end
    end
    
    %% The rest
    
    % Central difference velocity calculation
    time_sec_cum = cumsum(dtime_sec);
    vel_cd = zeros(1, length(dist_cum) - 2);
    for j=2:length(dist_cum)-1
        vel_cd(j)=(dist_cum(j+1)-dist_cum(j-1))/(time_sec_cum(j+1)-time_sec_cum(j-1));
    end
    
    v_smooth=smoothdata(vel_cd, 'gaussian', 10);
    ddist_vel_cd = dist_cum(2:end)-dist_cum(1:end-1);
    
    % Rolling Resistance Force and Energy
    F_r = Crr(idx) * m * g .* cos(theta);  % Newtons
    E_r = F_r .* ddist_xyz;  % Joules
    Er_cum = cumsum(E_r);  % Joules
    
    % Air density calculation
    rho_0 = 1.225;  % air density at sea level (kg/m^3)
    rho_1000 = 1.112;  % air density at 1000 m
    
    % air density interpolated varying with elevation
    rho = zeros(1, length(elevation));
    for i = 1:length(elevation)
        rho(i) = rho_0 + (rho_0-rho_1000)/(0 - 1000)*elevation(i);
    end
    
    % show density and elevation relationship
    figure(idx)
    subplot(2,1,1)
    plot(time_sec_cum, rho(1:end-1))
    xlabel('Time (s)')
    ylabel('Air Density (kg/m^3)')
    title('Air Density Variation')
    subplot(2,1,2)
    plot(time_sec_cum,elevation(1:end-1))
    xlabel('Time (s)')
    ylabel('Elevation (m)')
    title('Elevation Variation')
    
    linearEtravel_cum = zeros(1, len);
    linearE_tot = zeros(1, len);
    lagrangeEtravel_cum = zeros(1, len);
    lagrangeE_tot = zeros(1, len);
    
    for i=1:length(areaRange)
        linearA = areaRange(i);
        linearC_drag = linearDragVals(i);
        linearE_travel = 1/2 * linearC_drag * rho(i) * linearA * ddist_vel_cd .* v_smooth.^2;
        linearEtravel_cum(i) = sum(linearE_travel);
    
        linearE_tot(i) = Er_cum(end) + linearEtravel_cum(i);
    end

    for i=1:length(areaRange)
        lagrangeA = areaRange(i);
        lagrangeC_drag = lagrangeDragVals(i);
        lagrangeE_travel = 1/2 * lagrangeC_drag * rho(i) * lagrangeA * ddist_vel_cd .* v_smooth.^2;
        lagrangeEtravel_cum(i) = sum(lagrangeE_travel);
    
        lagrangeE_tot(i) = Er_cum(end) + lagrangeEtravel_cum(i);
    end
    
    linearE_tot = linearE_tot/1000;  % Joules -> kJ

    lagrangeE_tot = lagrangeE_tot/1000;  % J -> kJ
    
    % find index of the estimation thats closest to Strava's estimate
    [~, lineary] = min(abs(linearE_tot - stravaEst));
    
    % use index to get values
    linearAreaEst(idx) = areaRange(lineary);
    linearCDragEst(idx) = linearDragVals(lineary);
    linearEnergyExpenditureEst(idx) = linearE_tot(lineary);
    linearMPGEst(idx) = (Gas_yield/linearEnergyExpenditureEst(idx)) * dist_traveled ...
        * humanEnergyEff * milesPerMeter;
    linearPercentErr(idx) = (abs(stravaEst- linearEnergyExpenditureEst(idx))/stravaEst)*100;

    % find index of the estimation thats closest to Strava's estimate
    [~, lagrangey] = min(abs(linearE_tot - stravaEst));
    
    % use index to get values
    lagrangeAreaEst(idx) = areaRange(lagrangey);
    lagrangeCDragEst(idx) = lagrangeDragVals(lagrangey);
    lagrangeEnergyExpenditureEst(idx) = linearE_tot(lagrangey);
    lagrangeMPGEst(idx) = (Gas_yield/lagrangeEnergyExpenditureEst(idx)) * dist_traveled ...
        * humanEnergyEff * milesPerMeter;
    lagrangePercentErr(idx) = (abs(stravaEst- lagrangeEnergyExpenditureEst(idx))/stravaEst)*100;
end

% display results to command winder
disp("Linear Interpolation")
for idx=1:nRides
    ride = rides(idx);
    fprintf(['%s\nEstimated Frontal Area: %f m^2\nEstimated Energy Expenditure: %f kJ\n' ...
            'Estimated Energy Efficiency: %f mpg\nPercent Error: %f%%\n\n'], ...
            ride, linearAreaEst(idx), linearEnergyExpenditureEst(idx), linearMPGEst(idx), linearPercentErr(idx))
end

disp("Lagrange Interpolation")
for idx=1:nRides
    ride = rides(idx);
    fprintf(['%s\nEstimated Frontal Area: %f m^2\nEstimated Energy Expenditure: %f kJ\n' ...
            'Estimated Energy Efficiency: %f mpg\nPercent Error: %f%%\n\n'], ...
            ride, lagrangeAreaEst(idx), lagrangeEnergyExpenditureEst(idx), lagrangeMPGEst(idx), lagrangePercentErr(idx))
end

% Tabular results using lagrange interpolation
results = [rides'; lagrangeAreaEst; lagrangeCDragEst; ...
    lagrangeEnergyExpenditureEst; lagrangePercentErr; lagrangeMPGEst]';
results(3,:) = [];
resultsTable = array2table(results, "VariableNames", ...
    ["Ride", "Frontal Area (m^2)", "Drag Coefficient", ...
    "Energy Output (kJ)", "Relative Error to Stava (%)", "Efficiency (MPG)"]);
disp(resultsTable)
