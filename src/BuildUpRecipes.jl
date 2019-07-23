@userplot hornerplot

@recipe function f(j::hornerplot;
                    UnsteadyRange=false,Q=zeros(3), β=zeros(3), μ=ones(3), h=1,Rs=0.0, ct=1e-6, rw=0.583, ϕ=0.1, pr=true)

    T, P = j.args

    Tp=T[1]
    ΔT=T.-Tp
    HornerTime=(Tp.+ΔT)./ΔT
    H1hr=Tp.+1.0
    K=zeros(3)
        if UnsteadyRange!=false
       UnsteadyTime=HornerTime[(HornerTime.>=UnsteadyRange[1]) .& (HornerTime.<=UnsteadyRange[2])]
       UnsteadyPressure=P[(HornerTime.>=UnsteadyRange[1]) .& (HornerTime.<=UnsteadyRange[2])]

       m,b,r=Reg(UnsteadyTime,UnsteadyPressure,TypeReg="Logarithmic")
        UnsteadyLine=map(x->m*log10(x)+b,UnsteadyRange)
        UnsteadyExtrap=map(x->m*log10(x)+b,[UnsteadyRange[1],1])
       K[[1,2]]=map((q,be,mu)->(162.6*q*be*mu)/(abs(m)*h),Q[[1,2]],β[[1,2]],μ[[1,2]]) #liquis permeability Oil and Water
       K[3]= (162.6*(Q[3]*1000-Q[1]*Rs)*β[3]*μ[3])/(abs(m)*h)
       P1hr=m*log10(H1hr)+b
        λt=(K[1]/μ[1])+(K[2]/μ[2])+(K[3]/μ[3])
        S=1.151*(((P1hr-P[1])/abs(m))-(log10(λt/(ϕ*ct*rw^2)))+3.23)
        ΔPs=0.87*abs(m)*S
        Pi=b
        #Bg must be in bbl/scf
    end

    xlabel := "Horner Time"
    ylabel := "Pressure"
    xflip := true
    xscale := :log10
    legend := :topleft
    minorgrid := true
    gridcolor := :black
    ymirror := true

    @series  begin
        seriestype := :path
        linestyle := :dash
        linecolor := :black
        linewidth := 3
        label := "Horner Time"
        HornerTime, P
    end
    if  UnsteadyRange!=false
       @series begin
           seriestype := :path
           linestyle := :dash
           linecolor := :red
           linewidth := 1
           label := "Pi=$(round(Pi,digits=1)) psi"
           [UnsteadyRange[1],1],  UnsteadyExtrap
        end

        @series begin
            seriestype := :path
            linestyle := :solid
            linecolor := :red
            linewidth := 2
            label := "Unsteady Period"
            UnsteadyRange,  UnsteadyLine
         end
         pr==true ?  println("Unsteady State Parameters: \n Pi=$(round(Pi,digits=1)) psi  \n Ko=$(round(K[1],digits=1))md \n Kw=$(round(K[2],digits=1)) md \n Kg=$(round(K[3],digits=1)) md\n Skin=$(round(S,digits=2)) \n DeltaP Skin=$(round(ΔPs,digits=2)) psi ") : println(" ")
    end
end

@userplot derivativeplot

@recipe function f(j::derivativeplot)

        T, P = j.args

    Tp=T[1]
    Pwf=P[1]
    ΔT=T.-Tp
    ΔTe=ΔT./(1 .+(ΔT./Tp))
    ΔP=P.-Pwf
     δ=CentralDiff(ΔTe,ΔP)
    Der=((Tp.+ΔTe[2:end-1])./Tp).*δ.*ΔTe[2:end-1]


    xlabel := "Agarwal Time"
   ylabel := "Delta Pressure [psi]"
    legend := :topleft

    @series begin
        seriestype := :path
        linecolor := :black
        xscale := :log10
        yscale := :log10
        yformatter := :plain
        label := "Delta P"
        ΔTe[2:end], ΔP[2:end]
    end
      @series begin
        seriestype := :scatter
        markershape := :circle
        markersize := 2
        markercolor := :red
        xscale := :log10
        yscale := :log10
        yformatter := :plain
        label := "Derivative Function"
        ΔTe[2:end-1][Der.>0], Der[Der.>0]
    end

end
