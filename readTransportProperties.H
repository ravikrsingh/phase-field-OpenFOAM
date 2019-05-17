   IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

	dimensionedScalar pi = constant::mathematical::pi;
        dimensionedScalar dimt(transportProperties.lookup("dimt"));    //diemsion of time
        dimensionedScalar dimx(transportProperties.lookup("dimx"));    //diemsion of position
        dimensionedScalar m_1(transportProperties.lookup("m_1"));      //slope liquidus
        dimensionedScalar m_0(transportProperties.lookup("m_0"));      //slope solidus
        dimensionedScalar omega(transportProperties.lookup("omega"));  //relaxation coefficient for alpha (order parameter)
        dimensionedScalar gamma(transportProperties.lookup("gamma"));           //Surface Energy
        dimensionedScalar epsilon(transportProperties.lookup("epsilon"));       //Interface Width 
        dimensionedScalar c_Sol(transportProperties.lookup("c_Sol"));           //composition of solid in equilibrium with liquid
        dimensionedScalar c_Liq(transportProperties.lookup("c_Liq"));          //composition of liquid in quilibrium with solid
        dimensionedScalar c_eq(transportProperties.lookup("c_eq"));            //equilibrium composition or avg composition of alloy
        dimensionedScalar anti_trap(transportProperties.lookup("anti_trap"));           //antitrapping coefficient
        dimensionedScalar diff_Sol(transportProperties.lookup("diff_Sol"));             //diffusivity in solid
        dimensionedScalar diff_Liq(transportProperties.lookup("diff_Liq"));             //diffusivity in liquid
        dimensionedScalar G(transportProperties.lookup("G"));                           //thermal gradient
        dimensionedScalar v(transportProperties.lookup("v"));                           //velocity
        dimensionedScalar delta_01(transportProperties.lookup("delta_01"));             //strength of anisotropy
         dimensionedScalar A(transportProperties.lookup("A"));                          //G=Ac^2;d2G/dc^2=A
        dimensionedScalar T0(transportProperties.lookup("T0"));                         //Melting Temp
        dimensionedScalar noise_mag(transportProperties.lookup("noise_mag"));           //noise magnitude
        dimensionedScalar initial(transportProperties.lookup("initial"));               //constant value from temperature profile

