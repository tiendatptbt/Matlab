        a=[48.870398257397160;50;0.029039196408949;4.159676922543410e+03]
        P = a(1,1);
        I = a(2,1);
        D = a(3,1);
        sim('DC_SIM.slx');
        gravity = error'*error;