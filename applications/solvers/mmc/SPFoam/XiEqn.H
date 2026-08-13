dimensionedScalar tauRel("tauRel", dimTime, 0);

label deltaIndices = mesh.time().timeIndex() - mesh.time().startTimeIndex() - 1;

if(pdfMethodOn) // Added so solver can run without particle cloud
{
    tauRel.value() = max(pSets.first().coupling().tauRelaxStart().value() - deltaIndices*pSets.first().coupling().tauRelaxDelta(), pSets.first().coupling().tauRelaxTarget().value());

    if (pSets.first().coupling().tauUnits() == "timestep")
    {
        tauRel.value() = tauRel.value()*mesh.time().deltaT().value();
    }

    if (pSets.first().coupling().tauRelaxBlending())
    {
        Info << "tauRealaxStart = " << pSets.first().coupling().tauRelaxStart().value() << endl;
        Info << "tauRelaxTarget = " << pSets.first().coupling().tauRelaxTarget().value() << endl;
        Info << "tauRelaxDelta = " << pSets.first().coupling().tauRelaxDelta() << endl;
        Info << "tauRealax = " << tauRel.value() << endl;
    }
}

if (tauRel.value() < mesh.time().deltaT().value())
    tauRel.value() = mesh.time().deltaT().value();


tmp<fv::convectionScheme<scalar> > mvConvectionXi
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fieldsXi,
        phi,
        mesh.divScheme("div(phi,Xii)")
    )
);

forAll(Xi, XiI)
{
    volScalarField& Xii = Xi.Vars(XiI).field();
    volScalarField& XiiTarget = XiTarget.Vars(XiI).field();

    if (isA<couplingVar>(Xi.Vars(XiI)) || (Xi.Vars(XiI).refType()=="interpolated"))
    {
        // phiModEul (Eulerian phi-degree) is a particle-derived coupling field
        // with NO physical boundary conditions, so it cannot be a passively
        // transported scalar like z. Instead, relax it toward the particle-
        // projected phi-degree where the (second-conditioning) subset has
        // support, and leave the remaining cells with XiiTarget = Xii (no
        // relaxation) so they are filled by transport and diffusion. The target
        // is built once per time step (firstIter) and reused over PIMPLE.
        const bool isPhiModEul = (Xi.Vars(XiI).name() == "phiModEul");

        if (isPhiModEul)
        {
            if (pimple.firstIter())
            {
                XiiTarget = Xii;

                if (pdfMethodOn)
                {
                    scalarField& tgt = XiiTarget.primitiveFieldRef();

                    scalarField sumWt(mesh.nCells(), 0.0);
                    scalarField sumPhiWt(mesh.nCells(), 0.0);

                    const bool filterFlagged =
                        pSets.first().secondCondMixingEnabled();

                    // Same switch the KernelEstimation coupling reads: project
                    // the clean progress variable phi() when conditionOnPhi is
                    // set, else the OU-modulated phi-degree phiModified(). Read
                    // from the coupling coeffDict so the Eulerian phiModEul field
                    // and the kernel's particle coordinate stay consistent.
                    const bool condOnPhi =
                        pSets.first().coupling().coeffDict()
                            .lookupOrDefault<Switch>("conditionOnPhi", false);

                    forAllIters(pSets.first(), iter)
                    {
                        if (filterFlagged && iter().secondCondFlag() != 1)
                            continue;

                        const label celli = iter().cell();
                        const scalar w = iter().m();

                        sumWt[celli]    += w;
                        sumPhiWt[celli] +=
                            w * (condOnPhi ? iter().phi() : iter().phiModified());
                    }

                    forAll(tgt, celli)
                    {
                        if (sumWt[celli] > SMALL)
                            tgt[celli] = sumPhiWt[celli]/sumWt[celli];
                    }
                }
            }
        }
        else if (Xi.Vars(XiI).passive())
            XiiTarget = Xii;
        else
            XiiTarget = 0.0*Xii; // Need to add a way to calculate targets or
                                 // a source term

        fvScalarMatrix XiiEqn
        (
            fvm::ddt(rho,Xii)
//          + fvm::div(phi,Xii,"div(phi,Xii)")
          + mvConvectionXi->fvmDiv(phi, Xii)
          - fvm::laplacian(rho * DEff,Xii)
            ==
            rho * (XiiTarget - Xii) / tauRel
        );

        XiiEqn.relax();

        if (!pimple.finalIter())
        {
            XiiEqn.solve(mesh.solver("Xii"));
        }
        else
        {
            XiiEqn.solve(mesh.solver("XiiFinal"));
        }

        Xii.max(0.0);

        // phi-degree can exceed 1 (phi*exp(beta*omegaOU)); the upper clamp is
        // only valid for [0,1] coupling variables such as mixture fraction z.
        if (!isPhiModEul)
            Xii.min(1.0);
    }
    else
    {
        Info << "Internal Field is 0" <<endl;
        Xii.primitiveFieldRef() = 0.0;
    }
}
