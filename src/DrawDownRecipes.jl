@userplot drawdown

@recipe function f(j::drawdown;
                UnsteadyRange=false, PseudosteadyRange=false, pr=false,
                Q=zeros(3), β=zeros(3), μ=ones(3), h=1,Rs=0.0, ct=1e-6, rw=0.583, ϕ=0.1, B=1, To=0.1)

    T, P = j.args
    t=T[2:end]
    p=P[2:end]
    ΔP=P[1].-p

    Bsw=Q[2]/(Q[1]+Q[2])
    GOR=Q[3]*1000/Q[1]

    K=zeros(3)
    if UnsteadyRange!=false
       UnsteadyTime=t[(t.>=UnsteadyRange[1]) .& (t.<=UnsteadyRange[2])]
       UnsteadyPressure=p[(t.>=UnsteadyRange[1]) .& (t.<=UnsteadyRange[2])]

       m,b,r=Reg(UnsteadyTime,UnsteadyPressure,TypeReg="Logarithmic")
        UnsteadyLine=map(x->m*log10(x)+b,UnsteadyRange)
       K[[1,2]]=map((q,be,mu)->(162.6*q*be*mu)/(abs(m)*h),Q[[1,2]],β[[1,2]],μ[[1,2]]) #liquis permeability Oil and Water
       K[3]= (162.6*(Q[3]*1000-Q[1]*Rs)*β[3]*μ[3])/(abs(m)*h)
       P1hr=b
        λt=(K[1]/μ[1])+(K[2]/μ[2])+(K[3]/μ[3])
        S=1.151*(((P[1]-P1hr)/abs(m))-(log10(λt/(ϕ*ct*rw^2)))+3.23)
        ΔPs=0.87*abs(m)*S
        #Bg must be in bbl/scf
    end

    if PseudosteadyRange!=false
       PseudosteadyTime=t[(t.>=PseudosteadyRange[1]) .& (t.<=PseudosteadyRange[2])]
       PseudosteadyPressure=p[(t.>=PseudosteadyRange[1]) .& (t.<=PseudosteadyRange[2])]

       ms,bs,rs=Reg(PseudosteadyTime,PseudosteadyPressure,TypeReg="Linear")
        PseudosteadyLine=map(x->ms*x+bs,PseudosteadyRange)

        PV=(-sum(Q[1])*β[1])/(ms*ct*24*1e6)
        A=PV*1e6/(h*ϕ*7757.79)
        CA=5.456*(m/ms)*exp((2.303*(P1hr-bs))/m)
        rinv=0.0325*((K[1]*PseudosteadyRange[1])/(ϕ*μ[1]*ct))^0.5
    end

    #Log Log Plot Estimation Properties
    UnitLine=map(x->x*10^B,[t[1],To])

    C=(Q[1]*β[1]*To)/(24*UnitLine[2])
    WStime=((2e5+12000*S)*C)/(K[1]*h/μ[1])


    xlabel := "Time [hrs]"
    ylabel := "Pressure [psi]"
    xformatter := :auto
    layout := @layout [a;b c]
    size := (800,600)


    @series begin
       seriestype := :path
        linewidth := 3
       seriescolor := :darkgreen
        label := "DrawDown"
       subplot := 1
       t, p
    end

        if PseudosteadyRange!=false
        @series begin
        seriestype := :path
        linestyle := :solid
        linecolor := :red
        linewidth := 2
            subplot := 1
            label := "Pseudosteady Period"
      PseudosteadyRange, PseudosteadyLine
        end
    pr==true ?  println("Unsteady State Parameters: \n Pv=$(round(PV,digits=1)) MMbbl \n A=$(round(A,digits=1)) Acre\n CA=$(round(CA,digits=1)) \n Rinv=$(round(rinv,digits=1))  ") : println(" ")
    end
    @series begin
       seriestype := :path
        linestyle := :dashdotdot
       seriescolor := :darkgreen
       linewidth := 3
       xscale := :log10
        label := "DrawDown"
       subplot := 2
       t, p
    end
    if UnsteadyRange!=false
    @series begin
        seriestype := :path
        linestyle := :solid
        linecolor := :red
        linewidth := 2
            subplot := 2
            label := "Unsteady Period"
      UnsteadyRange, UnsteadyLine
        end

       pr==true ?  println("\nUnsteady State Parameters: \n Ko=$(round(K[1],digits=1)) md \n Kw=$(round(K[2],digits=1)) md \n Kg=$(round(K[3],digits=1)) md \n Skin=$(round(S,digits=2)) \n DeltaP Skin=$(round(ΔPs,digits=2)) psi ") : println(" ")
    end

    @series begin
        seriestype := :path
        linestyle := :solid
        linecolor := :black
        xscale := :log10
        yscale := :log10
        yformatter := :plain
        subplot :=3
        label := "Delta P psi"
        t,ΔP
    end

    @series begin
        seriestype := :path
        linestyle := :dash
        linecolor := :red
        minorgrid := true
        subplot := 3
        legend := :bottomright
        label := "Unit slope line"
        [t[1],To],UnitLine
    end
    pr==true ?  println("\nWellbore storage coefficient \n C=$(round(C,digits=2)) bbl/psi \nWellbore Storage End Time \n t=$(round(WStime,digits=2)) Hours") : println(" ")
end

@userplot ddsl

@recipe function f(j::ddsl;
        UnsteadyRange=false, pr=false,
                Q=zeros(3), β=zeros(3), μ=ones(3), h=1,Rs=0.0, ct=1e-6, rw=0.583, ϕ=0.1)

    T, P = j.args
    t=T[2:end]
    p=P[2:end]

    Bsw=Q[2]/(Q[1]+Q[2])
    GOR=Q[3]*1000/Q[1]

    K=zeros(3)
    if UnsteadyRange!=false
       UnsteadyTime=t[(t.>=UnsteadyRange[1]) .& (t.<=UnsteadyRange[2])]
       UnsteadyPressure=p[(t.>=UnsteadyRange[1]) .& (t.<=UnsteadyRange[2])]

       m,b,r=Reg(UnsteadyTime,UnsteadyPressure,TypeReg="Logarithmic")
        UnsteadyLine=map(x->m*log10(x)+b,UnsteadyRange)
       K[[1,2]]=map((q,be,mu)->(162.6*q*be*mu)/(abs(m)*h),Q[[1,2]],β[[1,2]],μ[[1,2]]) #liquis permeability Oil and Water
       K[3]= (162.6*(Q[3]*1000-Q[1]*Rs)*β[3]*μ[3])/(abs(m)*h)
       P1hr=b
        λt=(K[1]/μ[1])+(K[2]/μ[2])+(K[3]/μ[3])
        S=1.151*(((P[1]-P1hr)/abs(m))-(log10(λt/(ϕ*ct*rw^2)))+3.23)
        ΔPs=0.87*abs(m)*S
        #Bg must be in bbl/scf
    end

    xlabel := "Time [hrs]"
    ylabel := "Pressure [psi]"
    xformatter := :plain

    @series begin
       seriestype := :path
        linestyle := :dashdotdot
       seriescolor := :darkgreen
       linewidth := 3
       xscale := :log10
        label := "DrawDown"
       t, p
    end

    if UnsteadyRange!=false
    @series begin
        seriestype := :path
        linestyle := :solid
        linecolor := :red
        linewidth := 2

            label := "Unsteady Period"
      UnsteadyRange, UnsteadyLine
        end

       pr==true ?  println("Unsteady State Parameters: \n Ko=$(round(K[1],digits=1)) md \n Kw=$(round(K[2],digits=1)) md \n Kg=$(round(K[3],digits=1)) md \n Skin=$(round(S,digits=2)) \n DeltaP Skin=$(round(ΔPs,digits=2)) psi ") : println(" ")

    end

end

@userplot ddll

@recipe function f(j::ddll;
                    B=0.1, To=0.1, Q=0,β=1,  k=1, h=1, μ=1, s=0)

        T, P = j.args
    t=T[2:end]
    p=P[2:end]
    ΔP=P[1].-p

    UnitLine=map(x->x*10^B,[t[1],To])

    C=(Q*β*To)/(24*UnitLine[2])
    WStime=((2e5+12000*s)*C)/(k*h/μ)

    println("wellbore storage coefficient \n C=$(round(C,digits=2)) bbl/psi \nWellbore Storage End Time \n t=$(round(WStime,digits=2)) Hours")

    xlabel := "Time [hr]"
    ylabel := "Pressure [psi]"

    @series begin
        seriestype := :path
        linestyle := :solid
        linecolor := :black
        xscale := :log10
        yscale := :log10
        yformatter := :plain
        label := "Delta P psi"
        t,ΔP
    end

    @series begin
        seriestype := :path
        linestyle := :dash
        linecolor := :red
        minorgrid := true
        legend := :bottomright
        label := "Unit slope line"
        [t[1],To],UnitLine
    end
end

@userplot dd

@recipe function f(j::dd;
                    PseudosteadyRange=false, pr=false,
                   Q=zeros(3), β=zeros(3), μ=ones(3), h=1,Rs=0.0, ct=1e-6, rw=0.583, ϕ=0.1)
    T, P = j.args
    t=T[2:end]
    p=P[2:end]

        if PseudosteadyRange!=false
       PseudosteadyTime=t[(t.>=PseudosteadyRange[1]) .& (t.<=PseudosteadyRange[2])]
       PseudosteadyPressure=p[(t.>=PseudosteadyRange[1]) .& (t.<=PseudosteadyRange[2])]

       ms,bs,rs=Reg(PseudosteadyTime,PseudosteadyPressure,TypeReg="Linear")
        PseudosteadyLine=map(x->ms*x+bs,PseudosteadyRange)

        PV=(-sum(Q[1])*β[1])/(ms*ct*24*1e6)
        A=PV*1e6/(h*ϕ*7757.79)
    end

        xlabel := "Time [hr]"
    ylabel := "Pressure [psi]"

        @series begin
       seriestype := :path
        linewidth := 3
       seriescolor := :darkgreen
        label := "DrawDown"
       subplot := 1
       t, p
    end

        if PseudosteadyRange!=false
        @series begin
        seriestype := :path
        linestyle := :solid
        linecolor := :red
        linewidth := 2
            subplot := 1
            label := "Pseudosteady Period"
      PseudosteadyRange, PseudosteadyLine
        end
    println("Unsteady State Parameters: \n Pv=$(round(PV,digits=1)) MMbbl \n A=$(round(A,digits=1)) Acre")
    end
end
