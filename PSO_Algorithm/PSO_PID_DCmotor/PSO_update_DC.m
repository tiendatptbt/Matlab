
%Initial number
    %Init
% delta_x = 0.5;
% delta_y = 0.5;
% num_preg = 50;
clc
clear
clf

time=0;
num_g = 100;
num_p = 50;
    %Search range and Chose DIM of Planet
    % 3 DOF
    dim=3;
x_min = 0;
x_max = 100;
y_min = 0;
y_max = 100;
z_min = 0;
z_max = 10;
step = 10;
    %Create planet
x = (x_max-x_min)*rand(1,num_p)+x_min;
y = (y_max-y_min)*rand(1,num_p)+y_min;
z = (z_max-z_min)*rand(1,num_p)+z_min;
    %Init speed
w = 0.3;
% c1=[0.01 0;0 0.007];
% c2=[0.007 0;0 0.004];
c1 = 0.1;
c2 = 0.1;
    %Init hole
mini_hole = zeros(dim+1,1);
b_hole = zeros(dim+1,1);
    %Init velocity
v = zeros(dim,num_p);
val_min=0;
pos_min=0;
planet(1,:) = x;
planet(2,:) = y;
planet(3,:) = z;
% %Plot surface.
% [X1,X2] = meshgrid(x_min:step:x_max,y_min:step:y_max);
% C = -(X2+47).*sin(sqrt(abs(X2+X1/2+47)))-X1.*sin(sqrt(abs(X1-(X2+47))));
% %C = -(X2+47)*sin(sqrt(abs(X1/2+(X2+47))))-X1*sin(sqrt(abs(X1-(X2+47))));
% %C = -(X2+47)*sin(sqrt(abs(X1/2+(X2+47))))-X1*sin(sqrt(abs(X1-(X2+47))));
% %C = -abs(sin(X1).*cos(X2).*exp(abs(1-sqrt(X1.^2+X2.^2)/pi)));
% figure(2)
% subplot(3,1,1);
% surf(X1,X2,C);
% 
% subplot(3,1,2);
% contour(X1,X2,C,10);


% %Decade 2
for gen=1:num_g
    gen
    %Find the mini_hole
    for i=1:num_p
        P = planet(1,i);
        I = planet(2,i);
        D = planet(3,i);
%         P = 100;
%         I = 10;
%         D = 0;
        sim('REF_SIM.slx');
        gravity(1,i) = error'*error;%+0.3*u'*u;
    end
    [val_min,pos_min]=min_array(gravity);
    mini_hole(1:dim,1)=planet(:,pos_min);
    mini_hole(dim+1,1)=val_min;
    if gen==1
        b_hole(:,1) = mini_hole(:,1);
    end
    %Move each decade
    for j=1:num_p      
        v(:,j) = w*v(:,j) + c1*((b_hole(1:dim,1) - planet(:,j))) +c2*((mini_hole(1:dim,1) - planet(:,j)));
        planet(1,j)=planet(1,j)+v(1,j);%sqrt(v(1,j)^2+v(2,j)^2+v(3,j)^2)*cos(atan(v(2,j)/v(1,j))+(pi/1)*(0.5-rand(1)))*rand(1);%+v(1,j);
        planet(2,j)=planet(2,j)+v(2,j);%sqrt(v(1,j)^2+v(2,j)^2+v(3,j)^2)*sin(atan(v(2,j)/v(1,j))+(pi/1)*(0.5-rand(1)))*rand(1);%+v(2,j);
        planet(3,j)=planet(3,j)+v(3,j);%sqrt(v(1,j)^2+v(2,j)^2+v(3,j)^2)*sin(atan(v(2,j)/v(1,j))+(pi/1)*(0.5-rand(1)))*rand(1);%+v(3,j);
        % planet 1 limit
        if planet(1,j)>x_max
            planet(1,j)=x_max;
        elseif planet(1,j)<x_min
            planet(1,j)=x_min;
        end
        % planet 2 limit
        if planet(2,j)>y_max
            planet(2,j)=y_max;
        elseif planet(2,j)<y_min
            planet(2,j)=y_min;
        end
        % planet 3 limit
        if planet(3,j)>z_max
            planet(3,j)=z_max;
        elseif planet(3,j)<z_min
            planet(3,j)=z_min;
        end
            
    end
    %Find the b_hole each decade
    if (mini_hole(4,1)<=b_hole(4,1))
        b_hole(:,1) = mini_hole(:,1);
    end

% pause(0.2);

   
%     figure(2)
%     subplot(3,1,3);
%     plot(planet(1,:),planet(2,:),'.r')
%     axis([x_min x_max y_min y_max]);
%     subplot(3,1,2);
%     contour(X1,X2,C,20)
%     axis([x_min x_max y_min y_max]);
%     hold on
%     plot(mini_hole(1,1),mini_hole(2,1),'ob')
%     plot(b_hole(1,1),b_hole(2,1),'*r')
%     hold off
%     axis([x_min x_max y_min y_max]);
%     pause(0.2);
%     b_hole
end

figure(1)
subplot(2,1,1);
plot3(planet(1,:),planet(2,:),planet(3,:),'.r')
axis([0 50 0 100 0 0.1]);
subplot(2,1,2);
hold on
plot(gen,b_hole(4,1),'*r')
axis([0 10 0 500]);
hold off
P = b_hole(1);
I = b_hole(2);
D = b_hole(3);



