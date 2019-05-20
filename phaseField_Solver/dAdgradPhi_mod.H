     volVectorField q=dimx*fvc::grad(alpha);
     
     Info<< "Min/max magq:" << min(mag(q)).value() << ' '
     << max(mag(q)).value() << endl;
     

      ac_01 = 1.0 - delta_01*(3.0 - 4.0*(pow(q.component(0),4) + pow(q.component(1), 4) + pow(q.component(2), 4))/(1E-20 +pow(mag(q),4)));
     
     
      dAdq01 = 16.0*delta_01*((pow(q.component(0),3)/(1E-20 +pow(mag(q),4)) - q.component(0)*(pow(q.component(0),4) + pow(q.component(1), 4)+ pow(q.component(2), 4))/(1E-20 +pow(mag(q),6)))*vector(1,0,0)
      
      +(pow(q.component(1),3)/(1E-20 +pow(mag(q),4)) - q.component(1)*(pow(q.component(0),4) + pow(q.component(1), 4) + pow(q.component(2), 4))/(1E-20 +pow(mag(q),6)))*vector(0,1,0)
      
      +(pow(q.component(2),3)/(1E-20 +pow(mag(q),4)) - q.component(2)*(pow(q.component(0),4) + pow(q.component(1), 4) + pow(q.component(2), 4))/(1E-20 +pow(mag(q),6)))*vector(0,0,1));

 dadgradPhi = 2*gamma*(ac_01*dAdq01*mag(q)*mag(q));
 
 Info<< "Min/max ac_01:" << min(ac_01).value() << ' '
     << max(ac_01).value() << endl;
     
Info<< "Min/max dadgradPhi:" << min(dadgradPhi).value() << ' '
     << max(dadgradPhi).value() << endl;
