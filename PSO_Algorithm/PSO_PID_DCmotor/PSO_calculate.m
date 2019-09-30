function value= PSO_calculate(x,y)

    value = -(y+47)*sin(sqrt(abs(x/2+(y+47))))-x*sin(sqrt(abs(x-(y+47)))); %Eggholder function
    %value = sin(y)*exp((1-cos(x))^2)+cos(x)*exp((1-sin(y))^2)+(x-y)^2 ;   %Townsend function 
end