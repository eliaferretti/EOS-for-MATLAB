%
%  ______ ____   _____    _  _      __  __       _______ _               ____  
% |  ____/ __ \ / ____|  | || |    |  \/  |   /\|__   __| |        /\   |  _ \ 
% | |__ | |  | | (___    | || |_   | \  / |  /  \  | |  | |       /  \  | |_) |
% |  __|| |  | |\___ \   |__   _|  | |\/| | / /\ \ | |  | |      / /\ \ |  _ < 
% | |___| |__| |____) |     | |    | |  | |/ ____ \| |  | |____ / ____ \| |_) |
% |______\____/|_____/      |_|    |_|  |_/_/    \_\_|  |______/_/    \_\____/ 
%                                                                              
%                                                                              
%------------------------------------------------------------------------------------------
%
%
%	Politecnico di Milano (2022)
%
%	Bachelor's degree in chemical engineering
%
%	This code was developed and tested by Elia Ferretti
%
%	You can redistribute the code and/or modify it
%	Whenever the code is used to produce any publication or document,
%	reference to this work and author should be reported
%	No warranty of fitness for a particular purpose is offered
%	The user must assume the entire risk of using this code
%
%
%------------------------------------------------------------------------------------------




%Creator: ELIA FERRETTI Year: 2021 
 
%	Prerequisiti EoS 
%	Dichiarare globali le variabili: 
%		-  Tc = Vettore delle temperature critiche [K] 
%		-  pc = Vettore delle pressioni critiche [Pa] 
%		-  w = Vettore dei fattori acentrici di Pitzer [-] 
%
%	Prerequisiti gamma di Wilson 
%	Dichiarare globale la matrice dei coefficienti di Wilson con il nome “W” (doppia vu maiuscola)
% 
%	Significato dei termini di input 
%		-  T/Temp = Temperature [K] 
%		-  p/press = Pressione [Pa] 
%		-  x = composizione della fase (vettore delle frazioni molari) 
%		-  state = Stato fisico (“L” o “V”) 
%		-  i/index = indice del composto 
%
%	Note sui termini di output 
%		-  zeta [-] 
%		-  phi [-] 
%		-  entalpia residua [J/mol] 
%		-  gamma [-] 
%
%	Note 
%		-  Le funzioni sfruttano le unità di misura del Sistema Internazionale 
%		-  Tutte le funzioni sono state testate sui risultati numerici di esercizi proposti nel testo: 
%			“Fondamenti di Termodinamica dell’Ingegneria Chimica” – R.Rota 
 
 
 
  
function z = zetaVdW(index,T,p,state) 
global pc Tc w 
    R = 8.3144621; 
    a = 0.421875*R^2*Tc(index)^2/pc(index); 
    b = 0.125*R*Tc(index)/pc(index); 
   
    A = a*p/(R*T)^2; 
    B = b*p/(R*T); 
     
    polinomio = [1 -1-B A -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end    
end 
function z = zetaRK(index,T,p,state) 
global pc Tc w 
    R = 8.3144621; 
    k = (Tc(index)/T)^0.5; 
    a = 0.42748*R^2*Tc(index)^2*k/pc(index); 
    b = 0.08664*R*Tc(index)/pc(index); 
   
    A = a*p/(R*T)^2; 
    B = b*p/(R*T); 
     
    polinomio = [1 -1 A-B-B^2 -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end    
end 
function z = zetaRKS(index,T,p,state) 
global pc Tc w 
    R = 8.3144621; 
    S = 0.48+1.574*w(index)-0.176*w(index)^2; 
    k = (1+S(index)*(1-(T/Tc(index))^0.5))^2; 
    a = 0.42748*R^2*Tc(index)^2*k/pc(index); 
    b = 0.08664*R*Tc(index)/pc(index); 
   
    A = a*p/(R*T)^2; 
    B = b*p/(R*T); 
     
    polinomio = [1 -1 A-B-B^2 -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end    
end 
function z = zetaPR(index,T,p,state) 
global pc Tc w 
    R = 8.3144621; 
    S = 0.37464+1.54226*w(index)-0.26992*w(index)^2; 
    k = (1+S(index)*(1-(T/Tc(index))^0.5))^2; 
    a = 0.45724*R^2*Tc(index)^2*k/pc(index); 
    b = 0.0778*R*Tc(index)/pc(index); 
   
    A = a*p/(R*T)^2; 
    B = b*p/(R*T); 
     
    polinomio = [1 -1+B A-2*B-3*B^2 -A*B+B^2+B^3]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end    
end 
function z = zetaVIR(index,T,p) 
global pc Tc w 
    Tr = T/Tc(index); 
    pr = p/pc(index); 
    B0 = 0.083 - 0.422*Tr^(-1.6); 
    B1 = 0.139 - 0.172*Tr^(-4.2); 
    beta = B0 + B1*w(index); 
    z = 1 + beta*pr/Tr; 
end 
  
function hr = hrVdW(index,T,p,state) 
global Tc pc 
    R = 8.3144621; 
    a = 0.421875*R^2*Tc(index)^2/pc(index); 
    b = 0.125*R*Tc(index)/pc(index);     
    A = a*p/(R*T)^2; 
    B = b*p/(R*T); 
    polinomio = [1 -1-B A -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    hr = z-1-A/z; 
    hr = hr*R*T; 
end 
function hr = hrRK(index,T,p,state) 
global pc Tc w 
    R = 8.3144621; 
    k = (Tc(index)/T)^0.5; 
    a = 0.42748*R^2*Tc(index)^2*k/pc(index); 
    b = 0.08664*R*Tc(index)/pc(index); 
   
    A = a*p/(R*T)^2; 
    B = b*p/(R*T); 
     
    polinomio = [1 -1 A-B-B^2 -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    hr = z-1-3*A*log((z+B)/z)/(2*B); 
    hr = hr*R*T; 
end 
function hr = hrRKS(index,T,p,state) 
global pc Tc w 
    R = 8.3144621; 
    Tr = T/Tc(index); 
    S = 0.48+1.574*w(index)-0.176*w(index)^2; 
    k = (1+S*(1-(T/Tc(index))^0.5))^2; 
    a = 0.42748*R^2*Tc(index)^2*k/pc(index); 
    b = 0.08664*R*Tc(index)/pc(index); 
   
    A = a*p/(R*T)^2; 
    B = b*p/(R*T); 
     
    polinomio = [1 -1 A-B-B^2 -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    hr = z-1-A*(1+S*Tr^0.5/(1+S*(1-Tr^0.5)))*log((z+B)/z)/B; 
    hr = hr*R*T; 
end 
function hr = hrPR(index,T,p,state) 
global pc Tc w 
    R = 8.3144621; 
    Tr = T/Tc(index); 
    S = 0.37464+1.54226*w(index)-0.26992*w(index)^2; 
    k = (1+S*(1-(T/Tc(index))^0.5))^2; 
    a = 0.45724*R^2*Tc(index)^2*k/pc(index); 
    b = 0.0778*R*Tc(index)/pc(index); 
   
    A = a*p/(R*T)^2; 
    B = b*p/(R*T); 
     
    polinomio = [1 -1+B A-2*B-3*B^2 -A*B+B^2+B^3]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    hr = z-1-A*(1+S*Tr^0.5/(1+S*(1-Tr^0.5)))*log((z+B*(1+2^0.5))/(z+B*(1-2^0.5)))/(2^1.5*B); 
    hr = hr*R*T; 
end 
function hr = hrVIR(index,T,p) 
global pc Tc 
    R = 8.3144621; 
    Tr = T/Tc(index); 
    pr = p/pc(index); 
    B0 = 0.083 - 0.422*Tr^(-1.6); 
    B1 = 0.139 - 0.172*Tr^(-4.2); 
    beta = B0 + B1*w(index); 
    B = beta*R*Tc(index)/pc(index); 
    z = 1 + beta*pr/Tr; 
    dBdT = R*(0.675/Tr^2.6 + w(index)*0.722/Tr^5.2)/pc(index); 
    hr = p*(B-T*dBdT); 
end 
  
function phi = phiVdW(index,T,p,state) 
global Tc pc 
    R = 8.3144621; 
    a = 0.421875*R^2*Tc(index)^2/pc(index); 
    b = 0.125*R*Tc(index)/pc(index);     
    A = a*p/(R*T)^2; 
    B = b*p/(R*T); 
    polinomio = [1 -1-B A -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    phi = z-1-A/z-log(z-B); 
    phi = exp(phi); 
end 
function phi = phiRK(index,T,p,state) 
global pc Tc w 
    R = 8.3144621; 
    k = (Tc(index)/T)^0.5; 
    a = 0.42748*R^2*Tc(index)^2*k/pc(index); 
    b = 0.08664*R*Tc(index)/pc(index); 
   
    A = a*p/(R*T)^2; 
    B = b*p/(R*T); 
     
    polinomio = [1 -1 A-B-B^2 -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    phi = z-1-A*log((z+B)/z)/B-log(z-B); 
    phi = exp(phi); 
end 
function phi = phiRKS(index,T,p,state) 
global pc Tc w 
    R = 8.3144621; 
    Tr = T/Tc(index); 
    S = 0.48+1.574*w(index)-0.176*w(index)^2; 
    k = (1+S*(1-(T/Tc(index))^0.5))^2; 
    a = 0.42748*R^2*Tc(index)^2*k/pc(index); 
    b = 0.08664*R*Tc(index)/pc(index); 
   
    A = a*p/(R*T)^2; 
    B = b*p/(R*T); 
     
    polinomio = [1 -1 A-B-B^2 -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    phi = z-1-A*log((z+B)/z)/B-log(z-B); 
    phi = exp(phi); 
end 
function phi = phiPR(index,T,p,state) 
global pc Tc w 
    R = 8.3144621; 
    Tr = T/Tc(index); 
    S = 0.37464+1.54226*w(index)-0.26992*w(index)^2; 
    k = (1+S*(1-(T/Tc(index))^0.5))^2; 
    a = 0.45724*R^2*Tc(index)^2*k/pc(index); 
    b = 0.0778*R*Tc(index)/pc(index); 
   
    A = a*p/(R*T)^2; 
    B = b*p/(R*T); 
     
    polinomio = [1 -1+B A-2*B-3*B^2 -A*B+B^2+B^3]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    phi = z-1-A*log((z+B*(1+2^0.5))/(z+B*(1-2^0.5)))/(2^1.5*B)-log(z-B); 
    phi = exp(phi); 
end 
function phi = phiVIR(index,T,p) 
global pc Tc w 
    R = 8.3144621; 
    Tr = T/Tc(index); 
    pr = p/pc(index); 
    B0 = 0.083 - 0.422*Tr^(-1.6); 
    B1 = 0.139 - 0.172*Tr^(-4.2); 
    beta = B0 + B1*w(index); 
    B = beta*R*Tc(index)/pc(index); 
    phi = p*B/(R*T); 
    phi = exp(phi); 
end 

function y = zetaVdWmix(T,p,x,state) 
global pc Tc w 
amix = 0;bmix = 0; 
    R = 8.3144621; 
    for i=1:length(x) 
        a(i) = 0.421875*R^2*Tc(i)^2/pc(i); 
        b(i) = 0.125*R*Tc(i)/pc(i); 
    end 
  
    for i=1:length(x) 
        for j=1:length(x) 
            amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5; 
            bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2; 
        end 
    end 
    A = amix*p/(R*T)^2; 
    B = bmix*p/(R*T); 
     
    polinomio = [1 -1-B A -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    y = z; 
end 
function y = zetaRKmix(T,p,x,state) 
global pc Tc w 
amix = 0;bmix = 0; 
    R = 8.3144621; 
    for i=1:length(x) 
        Tr(i) = T/Tc(i); 
        k(i) = 1/Tr(i)^0.5; 
        a(i) = 0.42748*R^2*Tc(i)^2*k(i)/pc(i); 
        b(i) = 0.08664*R*Tc(i)/pc(i); 
    end 
  
    for i=1:length(x) 
        for j=1:length(x) 
            amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5; 
            bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2; 
        end 
    end 
    A = amix*p/(R*T)^2; 
    B = bmix*p/(R*T); 
     
    polinomio = [1 -1 A-B-B^2 -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    y = z;  
end 
function y = zetaRKSmix(T,p,x,state) 
global pc Tc w 
amix = 0;bmix = 0; 
    R = 8.3144621; 
    for i=1:length(x) 
        S(i) = 0.48+1.574*w(i)-0.176*w(i)^2; 
        Tr(i) = T/Tc(i); 
        k(i) = (1+S(i)*(1-Tr(i)^0.5))^2; 
  
        a(i) = 0.42748*R^2*Tc(i)^2*k(i)/pc(i); 
        b(i) = 0.08664*R*Tc(i)/pc(i); 
    end 
  
    for i=1:length(x) 
        for j=1:length(x) 
            amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5; 
            bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2; 
        end 
    end 
    A = amix*p/(R*T)^2; 
    B = bmix*p/(R*T); 
     
    polinomio = [1 -1 A-B-B^2 -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    y = z;  
end 
function y = zetaPRmix(T,p,x,state) 
global pc Tc w 
amix = 0;bmix = 0; 
    R = 8.3144621; 
    for i=1:length(x) 
        S(i) = 0.37464+1.54226*w(i)-0.26992*w(i)^2; 
        Tr(i) = T/Tc(i); 
        k(i) = (1+S(i)*(1-Tr(i)^0.5))^2; 
  
        a(i) = 0.45724*R^2*Tc(i)^2*k(i)/pc(i); 
        b(i) = 0.0778*R*Tc(i)/pc(i); 
    end 
  
    for i=1:length(x) 
        for j=1:length(x) 
            amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5; 
            bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2; 
        end 
    end 
    A = amix*p/(R*T)^2; 
    B = bmix*p/(R*T); 
     
    polinomio = [1 -1+B A-2*B-3*B^2 -A*B+B^2+B^3]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    y = z;  
end 
function y = zetaVIRmix(Temp,press,x) 
global pc Tc zc w 
    R = 8.3144621; 
    vc = R*Tc.*zc./pc; 
%Z 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            Z(i,j) = (zc(i)+zc(j))/2; 
        end 
    end 
%V     
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            V(i,j) = ((vc(i)^(1/3)+vc(j)^(1/3))/2)^3; 
        end 
    end 
%k - OKKIO KE NON SEMPRE C'E'   
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            k(i,j) = 1-(V(i,i)*V(j,j))^(1/2)/V(i,j); 
        end 
    end   
%T 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            T(i,j) = (Tc(i)*Tc(j))^(1/2)*(1-k(i,j)); 
        end 
    end 
%P 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            P(i,j) = Z(i,j)*R*T(i,j)/V(i,j); 
        end 
    end 
%W 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            W(i,j) = (w(i)+w(j))/2; 
        end 
    end 
    TR = Temp./T; 
    pR = press./P; 
     
    B0 = 0.083 - 0.422*TR.^(-1.6); 
    B1 = 0.139 - 0.172*TR.^(-4.2); 
    beta1 = B0 + B1.*W; 
    beta1 = beta1.*pR./TR; 
    beta = 0; 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            beta = beta + x(i)*x(j)*beta1(i,j); 
        end 
    end 
    z = 1 + beta; 
    y = z;  
end   

function hr = hrVdWmix(T,p,x,state) 
global pc Tc w 
amix = 0;bmix = 0; 
    R = 8.3144621; 
    for i=1:length(x) 
        a(i) = 0.421875*R^2*Tc(i)^2/pc(i); 
        b(i) = 0.125*R*Tc(i)/pc(i); 
    end 
  
    for i=1:length(x) 
        for j=1:length(x) 
            amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5; 
            bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2; 
        end 
    end 
    A = amix*p/(R*T)^2; 
    B = bmix*p/(R*T); 
     
    polinomio = [1 -1-B A -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    hr = z-1-A/z; 
    hr = hr*R*T; 
end 
function hr = hrRKmix(T,p,x,state) 
global pc Tc w 
amix = 0;bmix = 0; 
    R = 8.3144621; 
    for i=1:length(x) 
        Tr(i) = T/Tc(i); 
        k(i) = 1/Tr(i)^0.5; 
        a(i) = 0.42748*R^2*Tc(i)^2*k(i)/pc(i); 
        b(i) = 0.08664*R*Tc(i)/pc(i); 
    end 
  
    for i=1:length(x) 
        for j=1:length(x) 
            amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5; 
            bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2; 
        end 
    end 
    A = amix*p/(R*T)^2; 
    B = bmix*p/(R*T); 
     
    polinomio = [1 -1 A-B-B^2 -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    hr = z-1-3*A*log((z+B)/z)/(2*B); 
    hr = hr*R*T; 
end 
function hr = hrRKSmix(T,p,x,state) 
global pc Tc w 
amix = 0;bmix = 0; 
    R = 8.3144621; 
    for i=1:length(x) 
        S(i) = 0.48+1.574*w(i)-0.176*w(i)^2; 
        Tr(i) = T/Tc(i); 
        k(i) = (1+S(i)*(1-Tr(i)^0.5))^2; 
  
        a(i) = 0.42748*R^2*Tc(i)^2*k(i)/pc(i); 
        b(i) = 0.08664*R*Tc(i)/pc(i); 
    end 
  
    for i=1:length(x) 
        for j=1:length(x) 
            amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5; 
            bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2; 
        end 
    end 
    A = amix*p/(R*T)^2; 
    B = bmix*p/(R*T); 
     
    polinomio = [1 -1 A-B-B^2 -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    coeff = 0; 
    for i=1:length(x) 
        coeff = coeff + x(i)*S(i)*(a(i)*Tr(i)/k(i))^0.5; 
    end 
    hr = z-1-A*(1+coeff/amix^0.5)*log((z+B)/z)/B; 
    hr = hr*R*T; 
end 
function hr = hrPRmix(T,p,x,state) 
global pc Tc w 
amix = 0;bmix = 0; 
    R = 8.3144621; 
    for i=1:length(x) 
        S(i) = 0.37464+1.54226*w(i)-0.26992*w(i)^2; 
        Tr(i) = T/Tc(i); 
        k(i) = (1+S(i)*(1-Tr(i)^0.5))^2; 
  
        a(i) = 0.45724*R^2*Tc(i)^2*k(i)/pc(i); 
        b(i) = 0.0778*R*Tc(i)/pc(i); 
    end 
  
    for i=1:length(x) 
        for j=1:length(x) 
            amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5; 
            bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2; 
        end 
    end 
    A = amix*p/(R*T)^2; 
    B = bmix*p/(R*T); 
     
    polinomio = [1 -1+B A-2*B-3*B^2 -A*B+B^2+B^3]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    coeff = 0; 
    for i=1:length(x) 
        coeff = coeff + x(i)*S(i)*(a(i)*Tr(i)/k(i))^0.5; 
    end 
    hr = z-1-A*(1+coeff/amix^0.5)*log((z+B*(1+2^0.5))/(z+B*(1-2^0.5)))/(2^1.5*B); 
    hr = hr*R*T; 
end 
function hr = hrVIRmix(Temp,press,x) 
global pc Tc zc w 
    R = 8.3144621; 
    vc = R*Tc.*zc./pc; 
%Z 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            Z(i,j) = (zc(i)+zc(j))/2; 
        end 
    end 
%V     
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            V(i,j) = ((vc(i)^(1/3)+vc(j)^(1/3))/2)^3; 
        end 
    end 
%k - OKKIO KE NON SEMPRE C'E'   
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            k(i,j) = 1-(V(i,i)*V(j,j))^(1/2)/V(i,j); 
        end 
    end   
%T 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            T(i,j) = (Tc(i)*Tc(j))^(1/2)*(1-k(i,j)); 
        end 
    end 
%P 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            P(i,j) = Z(i,j)*R*T(i,j)/V(i,j); 
        end 
    end 
%W 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            W(i,j) = (w(i)+w(j))/2; 
        end 
    end 
    TR = Temp./T; 
    pR = press./P; 
     
    B0 = 0.083 - 0.422*TR.^(-1.6); 
    B1 = 0.139 - 0.172*TR.^(-4.2); 
    beta1 = B0 + B1.*W; 
    beta1 = beta1.*pR./TR; 
    beta = 0; 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            beta = beta + x(i)*x(j)*beta1(i,j); 
        end 
    end 
    z = 1 + beta; 
    B2 = 0.675./TR.^2.6; 
    B3 = 0.722./TR.^5.2; 
    dBdtpR = pR.*(B2 + W.*B3); 
    dBdt = 0; 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            dBdt = dBdt + x(i)*x(j)*dBdtpR(i,j); 
        end 
    end 
    hr = (beta - dBdt)*R*Temp; 
end 
  
function phi = phiVdWmix(T,p,x,index,state) 
global pc Tc w 
amix = 0;bmix = 0; 
    R = 8.3144621; 
    for i=1:length(x) 
        a(i) = 0.421875*R^2*Tc(i)^2/pc(i); 
        b(i) = 0.125*R*Tc(i)/pc(i); 
    end 
  
    for i=1:length(x) 
        for j=1:length(x) 
            amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5; 
            bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2; 
        end 
    end 
    A = amix*p/(R*T)^2; 
    B = bmix*p/(R*T); 
    Ai = a(index)*p/(R*T)^2; 
    Bi = b(index)*p/(R*T); 
     
    polinomio = [1 -1-B A -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    phi = Bi/(z-B)-2*(Ai*A)^0.5/z-log(z-B); 
    phi = exp(phi); 
end 
function phi = phiRKmix(T,p,x,index,state) 
global pc Tc w 
amix = 0;bmix = 0; 
    R = 8.3144621; 
    for i=1:length(x) 
        Tr(i) = T/Tc(i); 
        k(i) = 1/Tr(i)^0.5; 
        a(i) = 0.42748*R^2*Tc(i)^2*k(i)/pc(i); 
        b(i) = 0.08664*R*Tc(i)/pc(i); 
    end 
  
    for i=1:length(x) 
        for j=1:length(x) 
            amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5; 
            bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2; 
        end 
    end 
    A = amix*p/(R*T)^2; 
    B = bmix*p/(R*T); 
    Ai = a(index)*p/(R*T)^2; 
    Bi = b(index)*p/(R*T); 
     
    polinomio = [1 -1 A-B-B^2 -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    phi = Bi*(z-1)/B+A*(Bi/B-2*(Ai/A)^0.5)*log((z+B)/z)/B-log(z-B); 
    phi = exp(phi); 
end 
function phi = phiRKSmix(T,p,x,index,state) 
global pc Tc w 
amix = 0;bmix = 0; 
    R = 8.3144621; 
    for i=1:length(x) 
        S(i) = 0.48+1.574*w(i)-0.176*w(i)^2; 
        Tr(i) = T/Tc(i); 
        k(i) = (1+S(i)*(1-Tr(i)^0.5))^2; 
  
        a(i) = 0.42748*R^2*Tc(i)^2*k(i)/pc(i); 
        b(i) = 0.08664*R*Tc(i)/pc(i); 
    end 
  
    for i=1:length(x) 
        for j=1:length(x) 
            amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5; 
            bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2; 
        end 
    end 
    A = amix*p/(R*T)^2; 
    B = bmix*p/(R*T); 
    Ai = a(index)*p/(R*T)^2; 
    Bi = b(index)*p/(R*T); 
     
    polinomio = [1 -1 A-B-B^2 -A*B]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    coeff = 0; 
    for i=1:length(x) 
        coeff = coeff + x(i)*S(i)*(a(i)*Tr(i)/k(i))^0.5; 
    end 
    phi = Bi*(z-1)/B+A*(Bi/B-2*(Ai/A)^0.5)*log((z+B)/z)/B-log(z-B); 
    phi = exp(phi);  
end 
function phi = phiPRmix(T,p,x,index,state) 
global pc Tc w 
amix = 0;bmix = 0; 
    R = 8.3144621; 
    for i=1:length(x) 
        S(i) = 0.37464+1.54226*w(i)-0.26992*w(i)^2; 
        Tr(i) = T/Tc(i); 
        k(i) = (1+S(i)*(1-Tr(i)^0.5))^2; 
  
        a(i) = 0.45724*R^2*Tc(i)^2*k(i)/pc(i); 
        b(i) = 0.0778*R*Tc(i)/pc(i); 
    end 
  
    for i=1:length(x) 
        for j=1:length(x) 
            amix = amix + x(i)*x(j)*(a(i)*a(j))^0.5; 
            bmix = bmix + x(i)*x(j)*(b(i)+b(j))/2; 
        end 
    end 
    A = amix*p/(R*T)^2; 
    B = bmix*p/(R*T); 
    Ai = a(index)*p/(R*T)^2; 
    Bi = b(index)*p/(R*T); 
     
    polinomio = [1 -1+B A-2*B-3*B^2 -A*B+B^2+B^3]; 
    z = roots(polinomio); %cerco le radici 
    z = z(imag(z)==0); %elimino le radici immaginarie 
    if state=='L' 
        z = min(z); %z liquido  
    else 
        z = max(z); %z vapore 
    end 
    phi = Bi*(z-1)/B+A*(Bi/B-2*(Ai/A)^0.5)*log((z+B*(1+2^0.5))/(z+B*(1-2^0.5)))/(2^1.5*B)-log(z-B); 
    phi = exp(phi); 
end 
%OCCHIO A Kij 
function phi = phiVIRmix(Temp,press,x,index) 
global pc Tc w 
    R = 8.3144621; 
    zc = zeros(1,length(Tc)); 
     
    for i=1:length(Tc) 
        zc(i) = 1 - 0.339 - 0.033*w(i); 
    end 
    vc = R*Tc.*zc./pc; 
%Z 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            Z(i,j) = (zc(i)+zc(j))/2; 
        end 
    end 
%V     
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            V(i,j) = ((vc(i)^(1/3)+vc(j)^(1/3))/2)^3; 
        end 
    end 
%k - OKKIO KE NON SEMPRE C'E'   
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            k(i,j) = 1-(V(i,i)*V(j,j))^(1/2)/V(i,j); 
        end 
    end   
%T 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            T(i,j) = (Tc(i)*Tc(j))^(1/2)*(1-k(i,j)); 
        end 
    end 
%P 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            P(i,j) = Z(i,j)*R*T(i,j)/V(i,j); 
        end 
    end 
%W 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            W(i,j) = (w(i)+w(j))/2; 
        end 
    end 
    TR = Temp./T; 
    pR = press./P; 
     
    B0 = 0.083 - 0.422*TR.^(-1.6); 
    B1 = 0.139 - 0.172*TR.^(-4.2); 
    beta1 = B0 + B1.*W; 
    beta1 = beta1.*pR./TR; 
    beta = 0; 
    for i=1:length(Tc) 
        for j=1:length(Tc) 
            beta = beta + x(i)*x(j)*beta1(i,j); 
        end 
    end 
    beta2 = 0; 
    for i=1:length(Tc) 
        beta2 = beta2 + x(i)*beta1(index,i); 
    end 
    phi = -beta + 2*beta2; 
    phi = exp(phi); 
end 
  
%WILSON 
function y=gammaWilson(i,x) 
global W 
        beta = 0; 
        delta = 0; 
        for j=1:length(x) 
            beta = beta + x(j)*W(i,j); 
        end 
        for k=1:length(x) 
            eps = 0; 
            for j=1:length(x) 
                eps = eps + x(j)*W(k,j); 
            end 
            delta = delta + x(k)*W(k,i)/eps; 
        end 
        y = -log(beta)+1-delta; 
        y = exp(y); 
end 
 