clear all
R_Gaussian = 0.5:0.1:3;
E_Gaussian = [-1.0429963,-1.1011282,-1.117349,-1.1108504,-1.091914,-1.0661086,-1.0365389,-1.0051067,-0.9731106,-0.9414807, -0.9108736, -0.8817325, -0.8543376,-0.8288481,-0.8053328,-0.7837927,-0.7641777,-0.7464014,-0.7303533,-0.7159101,-0.7029436,-0.6913276,-0.6809408,-0.6716689,-0.6634047,-0.6560483];
error = 1e-12;
itr =0;
R=0.4:0.01:15;
for r = 0.4:0.01:15
    itr=itr+1;
    alpha1=[0.168856;0.6239133;3.42525];
    alpha2=[0.168856;0.6239133;3.42525];
    coeff1 =[0.444635;0.535328;0.154329];
    coeff2=[0.444635;0.535328;0.154329];
    ra =[0,r];
    S = zeros(2,2);
    for i = 1:2
        for j=1:2
    
            for k=1:3
                for l= 1:3 
                    p= alpha1(k)+alpha2(l);
                    q = alpha1(k)*alpha2(l)/p;
                    Rp = (alpha1(k)*ra(i)+alpha2(l)*ra(j))/(p);
                    S(i,j)= S(i,j)+ (p*q*4/pi/pi)^.75*(pi/p)^1.5*exp(-q*(ra(i)-ra(j))^2)*coeff1(k)*coeff2(l);
                    
                    
                    
                end
            end
           
        end
    end
    
    T=zeros(2,2);
    for i = 1:2
        for j=1:2
            
            for k=1:3
                for l= 1:3 
                    p= alpha1(k)+alpha2(l);
                    q = alpha1(k)*alpha2(l)/p;
                    Rp = (alpha1(k)*ra(i)+alpha2(l)*ra(j))/(p);
                    
                    
                    T(i,j)= T(i,j)+ (p*q*4/pi/pi)^.75*(pi/p)^1.5*q*(3-2*q*(ra(i)-ra(j))^2)*exp(-q*(ra(i)-ra(j))^2)*coeff1(k)*coeff2(l);
                    
                end
            end
           
        end
    end
    F0 = @(t) (pi/t)^0.5*erf(t^0.5)/2;
    Vne1=zeros(2,2);
    Vne2 = zeros(2,2);
    for i = 1:2
        for j=1:2
            
            for k=1:3
                for l= 1:3 
                    p= alpha1(k)+alpha2(l);
                    q = alpha1(k)*alpha2(l)/p;
                    Rp = (alpha1(k)*ra(i)+alpha2(l)*ra(j))/(p);
                    Rc1 = 0;
                    Rc2 = r;
                    
                    if (Rp-Rc1~=0)
                        Vne1(i,j)= Vne1(i,j)-2*pi* (p*q*4/pi/pi)^.75/p*F0(p*(Rp-Rc1)^2)*exp(-q*(ra(i)-ra(j))^2)*coeff1(k)*coeff2(l);
                    elseif (Rp==Rc1) 
                        Vne1(i,j)= Vne1(i,j)-2*pi* (p*q*4/pi/pi)^.75/p*exp(-q*(ra(i)-ra(j))^2)*coeff1(k)*coeff2(l);
                    end
                    if (Rp-Rc2~=0)
                    Vne2(i,j)= Vne2(i,j)-2*pi* (p*q*4/pi/pi)^.75/p*F0(p*(Rp-Rc2)^2)*exp(-q*(ra(i)-ra(j))^2)*coeff1(k)*coeff2(l);
                    elseif (Rp==Rc2)
                        Vne2(i,j)= Vne2(i,j)-2*pi* (p*q*4/pi/pi)^.75/p*exp(-q*(ra(i)-ra(j))^2)*coeff1(k)*coeff2(l);
                    end
    
    
                end
            end
           
        end
    end
    
    Vee=zeros(2,2,2,2);
    for i=1:2
        for j=1:2
            for k = 1:2
                for l=1:2
                    
                    for a=1:3
                        for b= 1:3 
                            for c=1:3
                                for d=1:3
                                    x= alpha1(a)+alpha2(b);
                                    y= alpha1(c)+alpha2(d);
    
                                    qx = alpha1(a)*alpha2(b)/x;
                                    qy= alpha1(c)*alpha2(d)/y;
                                    Rp = (alpha1(a)*ra(i)+alpha2(b)*ra(j))/(x);
                                    Rq = (alpha1(c)*ra(k)+alpha2(d)*ra(l))/(y);
                                    
                                    if (Rp-Rq~=0)
                                        Vee(i,j,k,l)= Vee(i,j,k,l)+2*pi^2.5/(x*y*(x+y)^0.5)*F0(x*y/(x+y)*(Rp-Rq)^2)*exp(-qx*(ra(i)-ra(j))^2-qy*(ra(k)-ra(l))^2)*coeff1(a)*coeff2(b)*coeff1(c)*coeff2(d)*(x*qx*y*qy*16/pi^4)^0.75;
                                    elseif (Rp==Rq) 
                                        Vee(i,j,k,l)= Vee(i,j,k,l)+2*pi^2.5/(x*y*(x+y)^0.5)*exp(-qx*(ra(i)-ra(j))^2-qy*(ra(k)-ra(l))^2)*coeff1(a)*coeff2(b)*coeff1(c)*coeff2(d)*(x*qx*y*qy*16/pi^4)^0.75;
                                    end
                                    
            
            
                               end
                            end
                   
                        end  
                    end
                end
            end
        end
    end
    
    Vnn= 1/(r);
    Hcore= T+Vne1+Vne2;
    Elast =1 ;
    count =0;
    
    C = ones(2,1);
    
    X=S^(-1/2);
    
    Eval(1)=0;
    
    while(abs(Elast - Eval(1))>error)
        P= 2.*C*C';
        count = count +1;
        G=zeros(2,2);
        Elast = Eval(1);
        for i=1:2
            for j=1:2
                for k=1:2
                    for l=1:2
                        G(i,j)=G(i,j)+P(k,l)*(Vee(i,j,l,k)-Vee(i,k,l,j)/2);
                    end
                end
            end
        end

        
        F=Hcore + G;
        Ft = X'*F*X;
        [V,val]=eig(Ft);
        eigval=diag(val,0);
        [Eval,index]=sort(eigval);
        Egr(count) = Eval(1);
        Cgr=X'*V(:,index(1));
        C = Cgr;
       
    
    end
    
    Ef =0;
    for i = 1:2
        for j = 1:2
            Ef =Ef + 0.5.*P(j,i)*(Hcore(i,j)+F(i,j));
        end
    end

    Ef =Ef + Vnn;
    E(itr)=Ef;
    
    if E(itr)==min(E)
        rmin = r;
        
        
        
    end
   

end

figure();
plot(R*0.529177,E,'b--','LineWidth',2);
title('Plot showing total energy from Hartree Fock  H-H  (From code vs Gaussian)');
xlabel('H-H bondlength (Angstrom)');
ylabel('Total energy (in a.u.)');

disp('The least energy bond length (A) for H-H is '); disp(rmin*0.529177);
disp('The energy (Hartree) of this bond is '); disp(min(E));

hold on;
plot(R_Gaussian,E_Gaussian,'m--','Linewidth',2);
legend('From code','From Gaussian');













