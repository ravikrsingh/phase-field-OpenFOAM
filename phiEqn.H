Random obj(1);
const scalar randNumber(obj.scalar01());
#include "dAdgradPhi_mod.H"
fvScalarMatrix alphaEqn
(
	omega*epsilon*dimt*fvm::ddt(alpha)     //alpha=1 liquid phase
                ==
		2.0*epsilon*gamma*dimx*dimx*fvm::laplacian(ac_01*ac_01, alpha) + epsilon*dimx*fvc::div(dadgradPhi)
                - 18.0*(gamma/epsilon)*(alpha)*(1-alpha)*(1-2.0*alpha)
                + 6*noise_mag*alpha*(1.0-alpha)*alpha*(1-alpha)*randNumber
                -(c_Sol-c_Liq)*(mu - (2*A*c_eq - 2*A*(T0-T)/m_1))*30.0*alpha*alpha*(1.0-alpha)*(1.0-alpha)

);

alphaEqn.solve();

volVectorField n=dimx*fvc::grad(alpha)/(1E-20+mag(dimx*fvc::grad(alpha)));
fvScalarMatrix muEqn
(
         (0.5)*dimt*fvm::ddt(mu) == diff_Liq*0.5*dimx*dimx*fvm::laplacian(alpha,mu) + (c_Sol-c_Liq)* dimt*fvc::ddt(alpha)*30.0*alpha*alpha*(1.0-alpha)*(1.0-alpha)- anti_trap*(c_Sol-c_Liq)*dimx*fvc::div((n*dimt*fvc::ddt(alpha)))
);

muEqn.solve();


    Info<< "Min/max alpha:" << min(alpha).value() << ' '
    << max(alpha).value() << endl;

