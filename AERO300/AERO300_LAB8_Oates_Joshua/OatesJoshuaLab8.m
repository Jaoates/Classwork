% Joshua Oates - lab 8

%% section 0 - clean up
clear all;
close all;
clc

%% section 1 - freq at mode 1
mode = 1;
n_pts = 100;
plots = "no";
Tf = .1;

freq1 = beamvibe(mode,n_pts,plots,Tf);

%% section 2 - freq at mode 2
mode = 1;
n_pts = 100;
plots = "no";
Tf = .1;

freq2 = beamvibe(mode,n_pts,plots,Tf);

%% section 3 - Positive Q
mode = 1;
n_pts = 100;
plots = "yes";
Tf = .001;
animationSpeed = 1;
Q = .002;

beamvibe(mode,n_pts,plots,Tf,animationSpeed,Q);

%% section 4 - freq at mode 2
mode = 2;
n_pts = 100;
plots = "no";
Tf = .1;

freq1 = beamvibe(mode,n_pts,plots,Tf);

%% section 5 - freq vs mode
clc
n_pts = 500;
plots = "no";
Tf = .1;
freq=[];
for i = 1:20
    freq = [freq,beamvibe(i,n_pts,plots,Tf)];
    disp(i)
end

%% section 6 - plot freq vs mode
figure
mode = 1:20;
semilogy(mode,freq)
legend("frequency","location","best")
xlabel("mode")
ylabel("frequency (Hz)")

