function [p] = fcn_ik_q_2_p(q)
% inverse kinematics: joint angles (q) --> actuator angles (p)
% q, p are 5x1 vectors of angles

p = zeros(5,1);

  p(1,1)=q(1);
  p(2,1)=atan2((8*cos(q(3))*(((7*sin(q(2)))/625 + (8*cos(q(2))*sin(q(3)))/625 - 57/2500)^2 + (64*...
         cos(q(3))^2)/390625 - ((49*cos(q(2)))/5000 + (399*sin(q(2)))/20000 + (57*cos(q(2))*sin(q(3)))/2500 - (7*sin(q(2))*...
         sin(q(3)))/625 - 3021/160000)^2)^(1/2))/625 - ((7*sin(q(2)))/625 + (8*cos(q(2))*sin(q(3)))/625 - 57/2500)*((49*cos(q(2)))/5000 + (399*...
         sin(q(2)))/20000 + (57*cos(q(2))*sin(q(3)))/2500 - (7*sin(q(2))*sin(q(3)))/625 - 3021/160000), - (8*cos(q(3))*((49*...
         cos(q(2)))/5000 + (399*sin(q(2)))/20000 + (57*cos(q(2))*sin(q(3)))/2500 - (7*sin(q(2))*...
         sin(q(3)))/625 - 3021/160000))/625 - ((7*sin(q(2)))/625 + (8*cos(q(2))*sin(q(3)))/625 - 57/2500)*(((7*sin(q(2)))/625 + (8*cos(q(2))*...
         sin(q(3)))/625 - 57/2500)^2 + (64*cos(q(3))^2)/390625 - ((49*cos(q(2)))/5000 + (399*sin(q(2)))/20000 + (57*cos(q(2))*...
         sin(q(3)))/2500 - (7*sin(q(2))*sin(q(3)))/625 - 3021/160000)^2)^(1/2));
  p(3,1)=atan2(((7*sin(q(2)))/625 - (8*cos(q(2))*sin(q(3)))/625 + 57/2500)*((49*...
         cos(q(2)))/5000 - (399*sin(q(2)))/20000 + (57*cos(q(2))*sin(q(3)))/2500 + (7*sin(q(2))*sin(q(3)))/625 - 3021/160000) + (8*...
         cos(q(3))*(((7*sin(q(2)))/625 - (8*cos(q(2))*sin(q(3)))/625 + 57/2500)^2 + (64*cos(q(3))^2)/390625 - ((49*...
         cos(q(2)))/5000 - (399*sin(q(2)))/20000 + (57*cos(q(2))*sin(q(3)))/2500 + (7*sin(q(2))*...
         sin(q(3)))/625 - 3021/160000)^2)^(1/2))/625, ((7*sin(q(2)))/625 - (8*cos(q(2))*sin(q(3)))/625 + 57/2500)*(((7*sin(q(2)))/625 - (8*cos(q(2))*...
         sin(q(3)))/625 + 57/2500)^2 + (64*cos(q(3))^2)/390625 - ((49*cos(q(2)))/5000 - (399*sin(q(2)))/20000 + (57*...
         cos(q(2))*sin(q(3)))/2500 + (7*sin(q(2))*sin(q(3)))/625 - 3021/160000)^2)^(1/2) - (8*cos(q(3))*((49*...
         cos(q(2)))/5000 - (399*sin(q(2)))/20000 + (57*cos(q(2))*sin(q(3)))/2500 + (7*sin(q(2))*sin(q(3)))/625 - 3021/160000))/625);
  p(4,1)=q(4) + (1007*pi)/1500 - atan2(((7*sin((19*pi)/30 + q(5)))/2500 + (21*sin(pi/9))/6250)*...
         ((147*cos((19*pi)/30 + q(5))*cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 +...
          (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 + 163349/6250000) - (((7*cos((19*pi)/30 + q(5)))/2500 + (21*...
         cos(pi/9))/6250 - 13/625)^2 - ((147*cos((19*pi)/30 + q(5))*cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 +...
          (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 + 163349/6250000)^2 + ((7*sin((19*pi)/30 + q(5)))/2500 +...
          (21*sin(pi/9))/6250)^2)^(1/2)*((7*cos((19*pi)/30 + q(5)))/2500 + (21*cos(pi/9))/6250 - 13/625), ((7*...
         cos((19*pi)/30 + q(5)))/2500 + (21*cos(pi/9))/6250 - 13/625)*((147*cos((19*pi)/30 + q(5))*...
         cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 + (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 +...
          163349/6250000) + ((7*sin((19*pi)/30 + q(5)))/2500 + (21*sin(pi/9))/6250)*(((7*cos((19*pi)/30 + q(5)))/2500 + (21*...
         cos(pi/9))/6250 - 13/625)^2 - ((147*cos((19*pi)/30 + q(5))*cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 +...
          (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 + 163349/6250000)^2 + ((7*sin((19*pi)/30 + q(5)))/2500 +...
          (21*sin(pi/9))/6250)^2)^(1/2));
  p(5,1)=q(4) - (1007*pi)/1500 + atan2(((7*sin((19*pi)/30 + q(5)))/2500 + (21*sin(pi/9))/6250)*...
         ((147*cos((19*pi)/30 + q(5))*cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 +...
          (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 + 163349/6250000) - (((7*cos((19*pi)/30 + q(5)))/2500 + (21*...
         cos(pi/9))/6250 - 13/625)^2 - ((147*cos((19*pi)/30 + q(5))*cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 +...
          (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 + 163349/6250000)^2 + ((7*sin((19*pi)/30 + q(5)))/2500 +...
          (21*sin(pi/9))/6250)^2)^(1/2)*((7*cos((19*pi)/30 + q(5)))/2500 + (21*cos(pi/9))/6250 - 13/625), ((7*...
         cos((19*pi)/30 + q(5)))/2500 + (21*cos(pi/9))/6250 - 13/625)*((147*cos((19*pi)/30 + q(5))*...
         cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 + (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 +...
          163349/6250000) + ((7*sin((19*pi)/30 + q(5)))/2500 + (21*sin(pi/9))/6250)*(((7*cos((19*pi)/30 + q(5)))/2500 + (21*...
         cos(pi/9))/6250 - 13/625)^2 - ((147*cos((19*pi)/30 + q(5))*cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 +...
          (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 + 163349/6250000)^2 + ((7*sin((19*pi)/30 + q(5)))/2500 +...
          (21*sin(pi/9))/6250)^2)^(1/2));

 