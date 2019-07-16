function CentralDiff(x,y)

    m=size(x,1)

    δ=zeros(m-2)
    for i=2:m-1
       δ[i-1]=(y[i-1]*((x[i]-x[i+1])/((x[i-1]-x[i])*(x[i-1]-x[i+1]))))+
              (y[i]*((2*x[i]-x[i-1]-x[i+1])/((x[i]-x[i-1])*(x[i]-x[i+1]))))+
              (y[i+1]*((x[i]-x[i-1])/((x[i+1]-x[i-1])*(x[i+1]-x[i]))))
    end
    return δ
end
