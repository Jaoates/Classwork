%% HW4 - Monte Carlo - Joshua Oates
clear all
close all
clc

addpath("C:\joshFunctionsMatlab")

%% Set up errors

err(1).name = "Local Sideral Time";
err(1).value = .01;
err(1).unit = "s";

err(2).name = "Target Location Latitude";
err(2).value = 10^-4;
err(2).unit = "deg";

err(3).name = "Target Location longitdue";
err(3).value = 10^-4;
err(3).unit = "deg";

err(4).name = "Target Location Elevation";
err(4).value = 10/1000;
err(4).unit = "km";

err(5).name = "Spacecraft Location Knowledge";
err(5).value = 3/1000;
err(5).unit = "km";

err(6).name = "Spacecraft Velo (in/cross track)";
err(6).value = .002/1000;
err(6).unit = "km/s";

err(7).name = "Spacecraft Velo (radial)";
err(7).value = .007/1000;
err(7).unit = "km/s";

err(8).name = "Sensor mounting location";
err(8).value = .01/1000;
err(8).unit = "km";

err(9).name = "Sensor mounting angle";
err(9).value = 10^-6;
err(9).unit = "deg";
%% geo parameters

r0 = [
    5980.8297
    -1282.3184
    4125.8019]; %km

v0 = [1.540887
    7.186813
    0]; %km/s

T0lla = [
    35.3
    -120.8
    .2]; % deg and km

jd = 2458981.1666; % days

n = 10000; % num runs

T0 = lla2eci_421(T0lla(1),T0lla(2),T0lla(3),jd + err(1).value*randn*1.15741e-5);
T0 = T0';

Re = norm(T0);
mue = 3.986004418 * (10^5);

%% run for default (part1)
disp("part1")
[e_vec, e_norm, Te, T] = geo_monte_sim(r0,v0,T0lla,jd,err,n);

e_norm = e_norm*1000; % m
e_norm = sort(e_norm);
e_90 = e_norm(floor(n*.9));
x = 1:length(e_norm);
figure
plot(x,e_norm)
title("sorted e-norm part1")

figure
scatter3(e_vec(1,:),e_vec(2,:),e_vec(3,:))
title("e-vec part1")




%% run for different altitudes (part2)
disp("part2")
r0hat = r0/norm(r0);
v0hat = v0/norm(v0);

m = 9;
e_90 =[];
for i = 1:m

    r0x = r0hat*((8000/(m-1))*i+Re);
    v0x = v0hat*sqrt(mue/norm(r0x));
    
    [e_vec, e_norm, Te, T] = geo_monte_sim(r0x,v0x,T0lla,jd,err,n);
    
    e_norm = e_norm*1000; % m
    e_norm = sort(e_norm);
    e_90(i) = e_norm(floor(n*.9));

    x = 1:length(e_norm);

%     figure
%     plot(x,e_norm)
%     title("sorted e-norm part2 & i ="+string(i))
%     
%     figure
%     scatter3(e_vec(1,:),e_vec(2,:),e_vec(3,:))
%     title("e-vec part2 & i ="+string(i))

end
i = 1:m;
figure
plot((8000/(m-1))*i+Re,e_90)
title("error90% vs alt")
xlabel("alt [km]")
ylabel("error90% in [m]")

%% latitude (part3)
disp("part3")
m = 10;% latitudes
e_90 = [];
for i = 1:m
    T0llax = T0lla;
    T0llax(1) = T0lla(1)+(10/m)*(i-1);
    [e_vec, e_norm, Te, T] = geo_monte_sim(r0x,v0x,T0lla,jd,err,n);
    e_norm = e_norm*1000; % m
    e_norm = sort(e_norm);
    e_90(i) = e_norm(floor(n*.9));

end
i = 1:m;
figure
plot(T0lla(1)+(10/m)*(i-1),e_90)
title("error90% vs lat")
xlabel("lat [deg]")
ylabel("error90% in [m]")

%% one var at a time (part4)
disp("part4")
%%%%%%%
% v0 = [0;0;0]
%%%%%%%

e_90 = [];
for j = 1:length(err)
    errx = err;
    for i = 1:length(err)
        errx(i).value = 0;
    end
    errx(j) = err(j);

    [e_vec, e_norm, Te, T] = geo_monte_sim(r0,v0,T0lla,jd,errx,n);
    
    e_norm = e_norm*1000; % m
    e_norm = sort(e_norm);
    e_90(j) = e_norm(floor(n*.9));

    x = 1:length(e_norm);

%     figure
%     plot(x,e_norm)
%     title("sorted e-norm part4 & i ="+string(j))
%     
%     figure
%     scatter3(e_vec(1,:),e_vec(2,:),e_vec(3,:))
%     title("e-vec part4 & i ="+string(j))
end

leg = categorical([err.name],[err.name]);
bar(leg,e_90)
xlabel("Error Type")
ylabel("Error Amount [m]")

disp("The ground error contributed by each part to the total error in m is as follows:")
for i = 1:length(err)
    disp(err(i).name+":")
    disp(string(e_90(i))+" m")
end

% part 5
disp("part5")
e_90byParts = rms(e_90);
disp("The ground error as estimated by the RMS of all errors above in m is:")
e_90byParts



%% part 6
disp("part6")
% Imager projected area
% assume alt of 1000km
FOV = 10e-3; %deg
z = 1000; % km
% this area side length
l = 2*sind(.5*FOV)*z;
l = l*1000; % m
% set r0 to be in the same dir as T0 and V0 to be perpendicular to it
T0hat = T0/norm(T0);
r0mag = (Re+1000);
r0x =T0hat*r0mag;
v0hat = joshFindPerpVec(r0x);
v0x = v0hat*sqrt(mue/norm(r0x));

% put this scenario into the thingy


[e_vec, e_norm, Te, T] = geo_monte_sim(r0x,v0x,T0lla,jd,err,n);

e_norm = e_norm*1000; % m
e_norm = sort(e_norm);
e_90 = e_norm(floor(n*.9));

% geometery worked out on paper...

l_90 = l-2*e_90;
lt = sqrt(.5*l_90^2);
At = lt^2;

disp("The area of the maximum target which can be imaged with 90% confidence in m^2 is")
At
disp("And has side length in m")
lt



%% dependancies
fList = matlab.codetools.requiredFilesAndProducts("C:\AERO421\HW4\HW4");
disp(fList')


