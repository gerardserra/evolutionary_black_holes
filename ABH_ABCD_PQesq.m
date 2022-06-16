function Var=ABH_ABCD_PQesq (fm,N,hring,xl,damping,geo,conf,rvec,dampconf,percent,rl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the reflection coefficient of an Acoustic Black Hole (ABH)
% termination using transfer matrices that relate the acoustic pressure (P)
% and the acoustic velocity (Q)
%
% Inputs:
%   fm = sampling frequency
%   N = number of rings
%   hring = ring thickness
%   xl = distance from the end termination for the truncation (m)
%   damping = damping value, typically 0.05
%   geo = 'a1' --> used to change the geometry. Typically set to the case named 'a1'
%   conf = 'lin','qua','vec' --> linear ABH, quadratic ABH or with the profile introduced by a vector
%   rvec = vector to stablish the inner radius of the black hole (in tan per cert with respect to R) for each ring.
%   dampconf = 'c0','kz' --> indicates if the damping is introduced in c0 or in kz
%                            (it is typically set to 'c0')
%   percent = tant per cent of resonators with porous material (only if 'kz' is set in dammpconf)
%   rl = []
%
% Outputs:
%   Var.R = reflection coefficient
%   Var.f = frequency vector
%   Var.alphan = absorption coefficient
%   Var.Rana = reflection coefficient from analytical solution
%   Var.Zin = input impedance
%   Var.k0L = k0L=k0*L;
%   Var.k0Ll = k0Ll=k0*(L-l);
%   Var.k0 = wavenumber
%   Var.L = length of the ABH termination (m)
%   Var.l = truncation length (m)
%   Var.c0 = speed of sound (m/s)
%   Var.N = number of tubes

% Examples for ABH optimization:
% >> ABH_ABCD_PQesq (1000,10,0.1/10,1.0e-3,0.05,'a1','vec',20:5:65,'kz',100); 

% Examples for linear and quad ABH:
% >> ABH_ABCD_PQesq (1000,40,0.001/1000,1.0e-3,0.05,'a1','lin',[],'kz',100);
% >> ABH_ABCD_PQesq (1000,40,0.5/100,1.0e-3,0.05,'a1','lin',[],'kz',100);
% >> ABH_ABCD_PQesq (1000,40,0.001/1000,1.0e-3,0.05,'a1','qua');
% >> ABH_ABCD_PQesq (1000,40,0.001/1000,1.0e-3,0.05,'a1','qua');
% >> ABH_ABCD_PQesq (1000,40,0.001/1000,1.0e-3,0.05,'a1','lin');
%
% Author: Marc Arnela, 2017
% Last update: 22 April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Codi               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global x1 x2 x3 x4

if exist('rl','var')==0
    rl=[];
end

if exist('dampconf','var')==0
    dampconf='c0';
end

if exist('rvec','var')==0
    rvec=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.- Geometria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('rl','var')~=0
    [L,R,hring,hcavity,Rr,xl,inici,l]=GeneraGeoABH(geo,N,hring,xl,rl,conf);
else
    [L,R,hring,hcavity,Rr,xl,inici,l]=GeneraGeoABH(geo,N,hring,xl,[],conf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.- Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0=340;
ro0=1.14;
rho0=ro0;
Z0=ro0*c0;

%fmax=fm/2;
fmax=2000;
NFFT=fm*50;

maxmostra = round(fmax * NFFT/fm);
factor = (maxmostra)/fmax;
f = [0:maxmostra]./factor;
f=f(1:end-1);
w=2*pi*f;
k=(w/c0)';
k0=k;

if exist('dampconf','var')==0
    c=c0.*(1+damping*1j);
    kz=(w./c).';
else
    if strcmp(dampconf,'kz')==0
        c=c0.*(1+damping*1j);
        kz=(w./c).';
    else
        mu=damping;
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3.- ACBD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=inici:N

%a) Paràmetres comuns
    [r1,r2,hcavity]=GeneraRadiusABH(conf,i,N,inici,xl,hcavity,hring,L,R,rl,rvec);
    %li=hcavity(i)+hring(i);
    
%b) Compensació per atenuació

    %kzcav --> propagation wave number in the cavity section
    %kzring --> propagation wave number in the ring section
    %kwcav --> wave number within the lumped element
    %Zzcav --> characteristic impedance in the propagation section of the cavity
    %Zzring --> characteristic impedance in the propagation section of the ring
    %Zwcav --> charactersitic impedance within the lumped element
    %Zl --> Wall impedance

    if strcmp(dampconf,'kz')==1
        if exist('percent','var')==0
            percent=100;
        end
        if (i<=percent/100*N)
            a3D=r2;
            if i==inici
                alphan(1:ceil(percent/100*N),1:length(f))=0;
            end
            [Zl,Zwcav,kwcav,alphan(i,:)]=CalculaImpedancia(f,rho0,c0,R-r2);
            mu=rho0*c0./Zl;
            %kzcav=k.*sqrt(1-1j.*(2*real(mu)./k/a3D)); %Sense aproximar l'arrel
            kzcav=k.*sqrt(1-1j.*(2*mu./k/a3D)); %Sense aproximar l'arrel
            Zzcav=Z0.*k./kzcav;
            kzring=k0;
            Zzring=Z0;
        else
            kzcav=k0;
            kzring=k0;
            kwcav=k0;
            Zzcav=Z0;
            Zzring=Z0;
            Zwcav=Z0;
        end
    else
        kzcav=kz;
        kzring=kz;
        kwcav=k0;
        Zzcav=Z0;
        Zzring=Z0;
        Zwcav=Z0;
    end   
    
%c) Coeficients Pressió/Cabal
    
    %Propagation
    m=pi*r2^2;
    [Apropring,Bpropring,Cpropring,Dpropring]=ABCD_ABH_Ring (hring(i),  kzring,Zzring,m); 
    [Apropcav,Bpropcav,Cpropcav,Dpropcav]    =ABCD_ABH_Ring (hcavity(i),kzcav, Zzcav,m);     
    [Aprop,Bprop,Cprop,Dprop]=ABCD_ABH_Product(Apropring,Bpropring,Cpropring,Dpropring,Apropcav,Bpropcav,Cpropcav,Dpropcav);   

    %Lumped
    [Alumpcav,Blumpcav,Clumpcav,Dlumpcav]=ABCD_ABH_Cavity (r1,r2,hcavity(i),R,kwcav,Zwcav);
   
    %Total
    [A2,B2,C2,D2]=ABCD_ABH_Product(Aprop,Bprop,Cprop,Dprop,Alumpcav,Blumpcav,Clumpcav,Dlumpcav);
    
%d) Càlcul matricial

    if(i==inici)
        Atot=A2;
        Btot=B2;
        Ctot=C2;
        Dtot=D2;
    else
        A1=Atot;
        B1=Btot;
        C1=Ctot;
        D1=Dtot;
                            
        [Atot,Btot,Ctot,Dtot]=ABCD_ABH_Product(A2,B2,C2,D2,A1,B1,C1,D1);                            
       
    end 
end

Zr=Inf;
Zin=(Atot + Btot ./ Zr)./(Ctot + Dtot./Zr);
ZinN=Zin/Z0*pi*R^2;
R=(ZinN-1)./(ZinN+1);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.- FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout==0
    %set(0,'DefaultTextInterpreter', 'tex');

    %fontsize = 14;
    %set(0,'defaultaxesfontsize',fontsize);
    %set(0,'defaulttextfontsize',fontsize);

    %figure;
    %subplot(211), plot(f/1000,abs(ZinN)); ylabel('$$|Z_{\rm in}|$$','Interpreter','latex');xlabel('Frequency (kHz)','Interpreter','latex'); ylim([0 10])
    %subplot(212), plot(f/1000,abs(R)); ylabel('$$|\mathcal{R}|$$','Interpreter','latex');  xlabel('Frequency (kHz)','Interpreter','latex'); ylim([0 1])
    
    %figure;
    %subplot(211), plot(k0*(L-l),abs(ZinN)); ylabel('$$|Z_{\rm in}|$$','Interpreter','latex');xlabh=xlabel('$$k_0 (L-l)$$','Interpreter','latex');set(xlabh,'Units','Normalized'); ylim([0 10]); xlim([0 4]); 
    %subplot(212), plot(k0*(L-l),abs(R)); ylabel('$$|\mathcal{R}|$$','Interpreter','latex');  xlabh=xlabel('$$k_0 (L-l)$$','Interpreter','latex');set(xlabh,'Units','Normalized'); ylim([0 1]);  xlim([0 4]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6.- OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(dampconf,'kz')==0
    switch conf(1:3)
        case 'lin'   
%             Rana=ABH_analytical_linear (k0,c',L,l);
        case 'qua'
%             Rana=ABH_analytical_quad   (k0,c',L,l);
    end
else
   Rana=[]; 
   c=[];
   Var.alphan=alphan;
end

k0Ll=k0*(L-l);
k0L=k0*L;

Var.f=f;
Var.R=R;
%Var.Rana=Rana;
Var.Zin=Zin;
Var.ZinN=ZinN;
Var.k0L=k0L;
Var.k0Ll=k0Ll;
Var.k0=k0;
Var.L=L;
Var.l=l;
Var.c0=c0;
Var.c=c;
Var.N=N;

end

function [Zl,Zw,kw,alphan]=CalculaImpedancia(f,rho0,c0,l)

paper='MechelBook'; %El que funciona millor

switch paper
    case 'Garai05'
        A=25.989;
        B=1.404;
        rhom=40; %Bulk density (kg/m^3)
%         l=40/1000;
        r=A*rhom^B; %airflow resistivity (Pa s/m^2)

        C1=0.078;
        C2=0.623;
        C3=0.074;
        C4=0.660;
        C5=0.159;
        C6=0.571;
        C7=0.121;
        C8=0.530;

        ZR=rho0*c0*          (1 + C1 * (rho0*f./r).^(-C2) );
        ZI=-rho0*c0 *        (C3 * (rho0*f/r).^(-C4) );
        alpha=(2*pi*f/c0) .* (C5 * (rho0*f/r).^(-C6) );
        beta= (2*pi*f/c0) .* (1 + C7 * (rho0*f./r).^(-C8) );
        
        Zw=ZR+1j*ZI;
        kw=alpha+1j*beta;    
        kw=kw.*(-1j); %Oriol Guasch after Scot1946
        
    case 'MechelBook'  %Porous media (based on Delany-Bazley 1970, and improved in Mechel 1976)
        
        r=1000;
        
        E=rho0*f/r; %Absorber variable E
        
        k0=2*pi*f/c0;
        
        Zw(1:length(E))=0;
        kw(1:length(E))=0;
        for i=1:length(E)
            if E(i)>1/60  %1/60!!!!
                Zw(i)=1+0.0485*E(i).^(-0.754)-1j*0.087*E(i).^(-0.73);
                Zw(i)=Zw(i)*rho0*c0;
                kw(i)=1-1j*0.189*E(i).^(-0.6185)+0.0978*E(i).^(-0.6929);
                kw(i)=kw(i).*k0(i);                
            else
                Zw(i)=(0.5/pi/E(i)+1j*1.4)  ./  sqrt(-1.466+1j*0.212/E(i))  ;
                Zw(i)=Zw(i)*rho0*c0;
                kw(i)=sqrt(1.466-1j*0.212./E(i));
                kw(i)=kw(i).*k0(i);           
            end
        end    
        
    case 'Munjal'  %Porous media pg 238 (based on Delany-Bazley 1970, and improved in Mechel 1976)
        
        r=1000;
        
        A=r./(rho0*f); %Normalitzed flow resistance of a lambda-deep layer
        
        k0=2*pi*f/c0;
        
        Zw(1:length(A))=0;
        kw(1:length(A))=0;
        for i=1:length(A)
            if A(i)<60
                Zw(i)=1+0.485.*A(i).^(0.754)-1j*0.087.*A(i).^(0.73);
                Zw(i)=Zw(i)*rho0*c0;
                kw(i)=-1j*0.189*A(i).^(0.6185)+1-1j*0.0978.*A(i).^(0.6929);
                kw(i)=kw(i).*k0(i);
            else
                Zw(i)=(0.5*A(i)+1j*1.4)./sqrt(-1.466+1j*0.212.*A(i));
                Zw(i)=Zw(i)*rho0*c0;
                kw(i)=sqrt(1.466-1j.*0.212*A(i));
                kw(i)=kw(i).*k0(i);           
            end
        end
        
    case 'Delany-Bazley_v1' %segons Garai05
        
        r=1000;
        
        C1=0.057;
        C2=0.754;
        C3=0.087;
        C4=0.732;
        C5=0.189;
        C6=0.595;
        C7=0.098;
        C8=0.700;

        ZR= rho0*c0*          (1 + C1 * (rho0*f./r).^(-C2) );
        ZI=-rho0*c0 *        (C3 * (rho0*f/r).^(-C4) );
        alpha=(2*pi*f/c0) .* (C5 * (rho0*f/r).^(-C6) );
        beta= (2*pi*f/c0) .* (1 + C7 * (rho0*f./r).^(-C8) );
        
        Zw=ZR+1j*ZI;
        kw=alpha+1j*beta;   
        kw=kw.*(-1j); %Oriol Guasch after Scot1946
       
        
    case 'Delany-Bazley_v2' %paper original 1970
        
        r=0.1;
        
        C1=9.08;
        C2=0.754;
        C3=11.9;
        C4=0.732;
        C5=10.3;
        C6=0.595;
        C7=10.8;
        C8=0.700;

        ZR= rho0*c0*          (1 + C1 * (f./r).^(-C2) );
        ZI=-rho0*c0 *        (C3 * (f/r).^(-C4) );
        alpha=(2*pi*f/c0) .* (C5 * (f/r).^(-C6) );
        beta= (2*pi*f/c0) .* (1 + C7 * (f./r).^(-C8) );
        
        Zw=ZR+1j*ZI;
        kw=alpha+1j*beta;       
        kw=kw.*(-1j); %Oriol Guasch after Scot1946
        
end

Zl=-1j*Zw.*cot(kw*l); %Munjal
ZlR=real(Zl);
ZlI=imag(Zl);

% alphan=(4*ZlR*rho0*c0)./(abs(Zl).^2 + 2*rho0*c0*ZlR + (rho0*c0).^2);

%alphan=1-abs((Zl-rho0*c0)./(Zl+rho0*c0)).^2;

 %hold all
 %figure(200);
 %subplot(211), plot(f,abs(Zl)); ylim([0 1000])
 %subplot(212), plot(f,abs(alphan));
 %plot(f,abs(alphan));

kw=kw.';
Zw=Zw.';
Zl=Zl.';

end

function [Atot,Btot,Ctot,Dtot]=ABCD_ABH_Product(A1,B1,C1,D1,A2,B2,C2,D2)
    Atot=A1.*A2+B1.*C2;
    Btot=A1.*B2+B1.*D2;
    Ctot=C1.*A2+D1.*C2;
    Dtot=C1.*B2+D1.*D2;  
end

function [A,B,C,D]=ABCD_ABH_Ring (li,kz,Z0,m)
    A=cos(kz*li);
    B=1j.*(Z0./m).*sin(kz*li);
    C=1j.*(m./Z0).*sin(kz*li);
    D=cos(kz*li); 
end

function [A,B,C,D]=ABCD_ABH_Cavity(r1,r2,hcavity,R,k0,Z0)

    Vcav=pi*hcavity*(R^2 - 1/3 * ( r2^2 + r1^2 + r2*r1  )  );
    Ycav=1j*k0*Vcav./Z0;
%     Scav=pi*(r2+r1)*sqrt(hcavity^2+(r2-r1)^2);      
    
%     S1=pi*r1*2;
%     S2=pi*r2*2;    
    
    A=1*ones(length(k0),1);
    B=0*ones(length(k0),1);
    C=Ycav;
    D=1*ones(length(k0),1);

%     A=1*ones(length(k0),1);
%     B=0*ones(length(k0),1);
%     C=0*ones(length(k0),1);
%     D=1*ones(length(k0),1);

end

function [L,R,hring,hcavity,Rr,xl,inici,l]=GeneraGeoABH (geo,N,hring,xl,rl,conf)

%a) Selecció geometria

switch geo
    
    case 'a1' %equiespaiat (ring constant), Rr=1
            
        L=0.5;
        R=0.23;        
                
        hring(1:N)=hring;        
        hcavity(1:N)=(L-N*hring)/N;
        inici=1;
        if exist('Rr','var')==0
            Rr=1;
        end           
        
    case 'a0' %equiespaiat (ring constant), Rr=0
            
        L=0.5;
        R=0.23;        
        
        hring(1:N)=hring;        
        hcavity(1:N)=(L-N*hring)/N;
        inici=1;
        if exist('Rr','var')==0
            Rr=0;
        end
                                                          
    case 'b1' %equiespaiat (cavity constant), Rr=1
        
        L=0.5;
        R=0.23;
        
        if exist('N','var')==0
            N=100;
        end

        hcavity(1:N)=hring;
        hring(1:N)=(L-N*hcavity(1))/N;

        if exist('Rr','var')==0
            Rr=1;
        end
        
    case 'c1' %Article internoise (variable)
        
        L=0.255;
        R=0.115;        
        N=18;
        xl=0;
        hring(1:N)=0.002;
        hcavity=ABCD_ABH_geometry_internoise(L,N,hring(1));
        inici=1;
        if exist('Rr','var')==0
            Rr=1;
        end
        
    case 'd1' %Article internoise (equiespaiat amb ring cavity constant)     
        
        L=0.255;
        R=0.115;        
        
        if exist('N','var')==0
            N=18;
        end

        xl=0;
        hring(1:N)=0.002;        
        hcavity(1:N)=(L-N*hring)/N;
        inici=1;
        if exist('Rr','var')==0
            Rr=1;
        end
        
    case 'e1' %Progressió lineal de cavities
        
        L=0.5;
        R=0.23;  
        
        if exist('N','var')==0
            N=8;
        end
            
        hring(1:N)=hring;
        hcavity=ABCD_ABH_geometry(L,N,R,hring(1));
        if exist('Rr','var')==0
            Rr=1;
        end     
        
    case 'f1' 
        
        L=0.5;
        R=0.23;  
        
        if exist('N','var')==0
            N=8;
        end
            
        hring(1:N)=hring;
        hcavity=ABCD_ABH_geometry2(L,N,hring(1));
        if exist('Rr','var')==0
            Rr=1;
        end             
      
end

%b) Control truncament

if length(conf)>4
    inici=1;
    xl=0;
    l=0;
else
    if isempty(rl)~=1
        switch conf(1:3)
            case 'lin'
                xl=rl/R*L;
            case 'qua'
                xl=sqrt(rl/R*L^2);
        end 
        display(['xl modifeid to ',num2str(xl)]);    
    end   

    if hcavity(1)==0 || hcavity(1)<0
        display('warning: hcavity<=0');
        Var=0;
        return
    end
    if (xl>hcavity(1))
        display('warning: xl>hcavity(1). ABH truncated.');
        aux=0;
        i=0;
        while aux<=xl
        %while (aux-xl)>hcavity(1)
            i=i+1;
            aux=aux+hcavity(i)+hring(i);    
        end
        inici=i;
        l=xl;
        xl=aux-xl;
        if xl>hcavity(inici)
            display('error: xl>hcavity(1) so the ring should be cut');
            hcavityI=hcavity(inici);        
            display(['inici=',num2str(inici),' xl=',num2str(xl),' hcavity=',num2str(hcavityI)]);  
            Var=0;
            return
        end
    else
        inici=1;
        l=xl;
    end
end
    
end

function [r1_out,r2_out,hcavity]=GeneraRadiusABH(conf,i,N,inici,xl,hcavity,hring,L,R,rl,rvec)

global x1 x2 x3 x4

    if i==inici
%         h100=figure(100);
%         close (h100);
%         h100=figure(100);
        xlim([-L 0]);
        if i==1
            x1(i)=-xl;
            x2(i)=-hcavity(i);        
            x3(i)=-hcavity(i)-hring(i)/2;                
            x4(i)=-hcavity(i)-hring(i);                        
        else
            x4(1)=-hcavity(1)-hring(1);                        
            for i2=2:inici
                x4(i2)=-hcavity(i2)-hring(i)+x4(i2-1);
            end  
            x1(inici)=x4(inici-1)-xl;
            x2(inici)=-hcavity(inici)+x4(inici-1);
            x3(inici)=-hcavity(inici)-hring(i)/2+x4(inici-1);                   
        end
        hcavity(i)=hcavity(i)-xl;
    else
        x1(i)=x4(i-1);
        x2(i)=-hcavity(i)+x4(i-1);
        x3(i)=-hcavity(i)-hring(i)/2+x4(i-1);
        x4(i)=-hcavity(i)-hring(i)+x4(i-1);
    end
           
    %b) Radius
    
    if length(conf)>3 %rlcorrected
        switch conf(1:3)
            case 'lin'
                if i==inici
                    r1(i)=-(R-rl)/L*x1(i)+rl;
                    r3(i)=-(R-rl)/L*x3(i)+rl;    
                    r2(i)=r3(i);    
                    r4(i)=r3(i);
                elseif i==N
                    r1(i)=-(R-rl)/L*(x1(i)+hring(i)/2)+rl;        
                    r2(i)=R;
                    r3(i)=R;        
                    r4(i)=R;                
                else
                    r1(i)=-(R-rl)/L*(x1(i)+hring(i)/2)+rl;
                    r3(i)=-(R-rl)/L*x3(i)+rl;    
                    r2(i)=r3(i);    
                    r4(i)=r3(i);        
                end    
            case 'qua'
                if i==inici
                    r1(i)=(R-rl)/L^2*x1(i)^2+rl;
                    r3(i)=(R-rl)/L^2*x3(i)^2+rl;    
                    r2(i)=r3(i);    
                    r4(i)=r3(i);
                elseif i==N
                    r1(i)=(R-rl)/L^2*(x1(i)+hring(i)/2)^2+rl;        
                    r2(i)=R;
                    r3(i)=R;        
                    r4(i)=R;                
                else
                    r1(i)=(R-rl)/L^2*(x1(i)+hring(i)/2)^2+rl;
                    r3(i)=(R-rl)/L^2*x3(i)^2+rl;    
                    r2(i)=r3(i);    
                    r4(i)=r3(i);        
                end    
        end
    else
        switch conf(1:3)
            case 'lin'   
                if i==inici
                    r1(i)=-R/L*x1(i);
                    r3(i)=-R/L*x3(i);    
                    r2(i)=r3(i);    
                    r4(i)=r3(i);
                elseif i==N
                    r1(i)=-R/L*(x1(i)+hring(i)/2);        
                    r2(i)=R;
                    r3(i)=R;        
                    r4(i)=R;                
                else
                    r1(i)=-R/L*(x1(i)+hring(i)/2);
                    r3(i)=-R/L*x3(i);    
                    r2(i)=r3(i);    
                    r4(i)=r3(i);        
                end
            case 'qua'            
                if i==inici
                    r1(i)=R/L^2*x1(i)^2;
                    r3(i)=R/L^2*x3(i)^2;    
                    r2(i)=r3(i);    
                    r4(i)=r3(i);
                elseif i==N
                    r1(i)=R/L^2*(x1(i)+hring(i)/2)^2;        
                    r2(i)=R;
                    r3(i)=R;        
                    r4(i)=R;                
                else
                    r1(i)=R/L^2*(x1(i)+hring(i)/2)^2;
                    r3(i)=R/L^2*x3(i)^2;    
                    r2(i)=r3(i);    
                    r4(i)=r3(i);        
                end
            case 'vec'   
                if i==inici
                    r1(i)=rvec(i)*R/100;
                    r3(i)=rvec(i+1)*R/100;    
                    r2(i)=r3(i);    
                    r4(i)=r3(i);
                elseif i==N
                    r1(i)=rvec(i)*R/100;        
                    r2(i)=R;
                    r3(i)=R;        
                    r4(i)=R;                
                else
                    r1(i)=rvec(i)*R/100;
                    r3(i)=rvec(i+1)*R/100;    
                    r2(i)=r3(i);    
                    r4(i)=r3(i);        
                end                
        end
    end
    
    
    % hold on
    % box on
    % xv=[x1(i) x1(i) x2(i) x2(i) x3(i) x4(i)];
    % yv=[r1(i) R     R     r2(i) r3(i) r4(i)];     
    % plot(xv,yv,'Color',	[0, 0.4470, 0.7410])
    % ylabel('ABH profile (m)','Interpreter','latex');
    % xlabel('Distance from termination (m)','Interpreter','latex');
    % hold off
    % drawnow

    r1_out=r1(i);
    r2_out=r2(i);
    
end