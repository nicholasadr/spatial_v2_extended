%TODO(@nic) had to convert model to struct:
m = struct(mit_humanoid_model());

N = 10;
t = linspace(0,1,N);
x = zeros(m.NQ,N);
x(1,:) = 1;
%TODO(@nicholasadr): these are link joints
%x(11,:) = linspace(0,2,N);
%x(12,:) = linspace(0,2,N);
showmotion(m,t,x);

%TODO(@nicholasadr): model.camera.direction in tello_model() is not working?
