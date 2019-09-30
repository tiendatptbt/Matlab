
%Initial number
    %Init
% delta_x = 0.5;
% delta_y = 0.5;
% num_preg = 50;
num_g = 100;
num_p = 100;
    %Search range
x_min = -512;
x_max = 512;
y_min = -512;
y_max = 512;
step = 10;
    %Create planet
x = (x_max-x_min)*rand(1,num_p)+x_min;
y = (y_max-y_min)*rand(1,num_p)+y_min;
    %Init speed
w = 0.6;
% c1=[0.01 0;0 0.007];
% c2=[0.007 0;0 0.004];
c1 = 0.1;
c2 = 0.1;
    %Init hole
mini_hole = zeros(3,1);
b_hole = zeros(3,1);
    %Init velocity
v = zeros(2,num_p);
val_min=0;
pos_min=0;
planet(1,:) = x;
planet(2,:) = y;
%Plot surface.
[X1,X2] = meshgrid(x_min:step:x_max,y_min:step:y_max);
C = -(X2+47).*sin(sqrt(abs(X2+X1/2+47)))-X1.*sin(sqrt(abs(X1-(X2+47))));
%C = -(X2+47)*sin(sqrt(abs(X1/2+(X2+47))))-X1*sin(sqrt(abs(X1-(X2+47))));
%C = -(X2+47)*sin(sqrt(abs(X1/2+(X2+47))))-X1*sin(sqrt(abs(X1-(X2+47))));
%C = -abs(sin(X1).*cos(X2).*exp(abs(1-sqrt(X1.^2+X2.^2)/pi)));
figure(2)
subplot(3,1,1);
surf(X1,X2,C);

subplot(3,1,2);
contour(X1,X2,C,10);

%Decade 1
% for pre_gen=1:num_preg
%     pre_gen
%     for i_pre=1:num_p
%         %s_1
%         satelite(1,1)=PSO_calculate(planet(1,i_pre),planet(2,i_pre));
%         satelite(2:3,1)=[planet(1,i_pre),planet(2,i_pre)];
%         %s_2
%         satelite(1,2)=PSO_calculate(planet(1,i_pre)+delta_x,planet(2,i_pre));
%         satelite(2:3,2)=[planet(1,i_pre)+delta_x,planet(2,i_pre)];
%         %s_3
%         satelite(1,3)=PSO_calculate(planet(1,i_pre)+delta_x,planet(2,i_pre)+delta_y);
%         satelite(2:3,3)=[planet(1,i_pre)+delta_x,planet(2,i_pre)+delta_y];
%         %s_4
%         satelite(1,4)=PSO_calculate(planet(1,i_pre),planet(2,i_pre)+delta_y);
%         satelite(2:3,4)=[planet(1,i_pre),planet(2,i_pre)+delta_y];
%         %s_5
%         satelite(1,5)=PSO_calculate(planet(1,i_pre)-delta_x,planet(2,i_pre)+delta_y);
%         satelite(2:3,5)=[planet(1,i_pre)-delta_x,planet(2,i_pre)+delta_y];
%         %s_6
%         satelite(1,6)=PSO_calculate(planet(1,i_pre)-delta_x,planet(2,i_pre));
%         satelite(2:3,6)=[planet(1,i_pre)-delta_x,planet(2,i_pre)];
%         %s_7
%         satelite(1,7)=PSO_calculate(planet(1,i_pre)-delta_x,planet(2,i_pre)-delta_y);
%         satelite(2:3,7)=[planet(1,i_pre)-delta_x,planet(2,i_pre)-delta_y];
%         %s_8
%         satelite(1,8)=PSO_calculate(planet(1,i_pre)-delta_x,planet(2,i_pre)-delta_y);
%         satelite(2:3,8)=[planet(1,i_pre)-delta_x,planet(2,i_pre)-delta_y];
%         %s_9
%         satelite(1,9)=PSO_calculate(planet(1,i_pre)+delta_x,planet(2,i_pre)-delta_y);
%         satelite(2:3,9)=[planet(1,i_pre)+delta_x,planet(2,i_pre)-delta_y];
%         p_move =1 ;
%         min_sa = satelite(1,1);
%         for si=2:9
%              if satelite(1,si) < min_sa
%                  p_move = si;
%                  min_sa = x(1,si);
%              end
%         end
%         planet(:,i_pre) = satelite(2:3,p_move);      
%     end
%     figure(1)
%     plot(planet(1,:),planet(2,:),'.r')
%      axis([-512 512 -512 512]);
% %    axis([-10 0 -6.5 0]);
%     pause(0.2);
% end
%Decade 2
for gen=1:num_g
    gen
    %Find the mini_hole
    for i=1:num_p
%         P = planet(1,i);
%         I = planet(2,i);
%         D = planet(3,i);
%         sim('DC_SIM.slx');
%         gravity(1,i) = e'*e;
       gravity(1,i)=PSO_calculate(planet(1,i),planet(2,i));
    end
    [val_min,pos_min]=min_array(gravity);
    mini_hole(1:2,1)=planet(:,pos_min);
    mini_hole(3,1)=val_min;
    if gen==1
        b_hole(:,1) = mini_hole(:,1);
    end
    %Move each decade
    for j=1:num_p      
        v(:,j) = w*v(:,j) + c1*((b_hole(1:2,1) - planet(:,j))) +c2*((mini_hole(1:2,1) - planet(:,j)));
        planet(1,j)=planet(1,j)+sqrt(v(1,j)^2+v(2,j)^2)*cos(atan(v(2,j)/v(1,j))+(pi/1)*(0.5-rand(1)))*rand(1);%+v(1,j);
        planet(2,j)=planet(2,j)+sqrt(v(1,j)^2+v(2,j)^2)*sin(atan(v(2,j)/v(1,j))+(pi/1)*(0.5-rand(1)))*rand(1);%+v(2,j);
        if planet(1,j)>x_max
            planet(1,j)=x_max;
        elseif planet(1,j)<x_min
            planet(1,j)=x_min;
        end
        if planet(2,j)>y_max
            planet(2,j)=y_max;
        elseif planet(2,j)<y_min
            planet(2,j)=y_min;
        end
            
    end
    %Find the b_hole each decade
    if (mini_hole(3,1)<=b_hole(3,1))
        b_hole(:,1) = mini_hole(:,1);
    end

   
    figure(2)
    subplot(3,1,3);
    plot(planet(1,:),planet(2,:),'.r')
    axis([x_min x_max y_min y_max]);
    subplot(3,1,2);
    contour(X1,X2,C,20)
    axis([x_min x_max y_min y_max]);
    hold on
    plot(mini_hole(1,1),mini_hole(2,1),'ob')
    plot(b_hole(1,1),b_hole(2,1),'*r')
    hold off
    axis([x_min x_max y_min y_max]);
    pause(0.2);
    b_hole
end



