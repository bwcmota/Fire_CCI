function [diamax valor valor_dist dj_max dj_min]=classifica_cor(sat, site, nir_cor, nir_ref, mask_gl, diajul);
intdias=clima(site);

alpha=0.05;
grau=1;
if sat==1
    numero=10;    
elseif sat==2   
    numero=6;
end

ref2=nir_ref';
ref=nir_cor';
mask=mask_gl';
validos=ref(find(mask~=1));
Ind=find(mask~=1);

diamax=0;
dj_max=0;
dj_min=0;
valor=0;
    
tam2=size(ref);
tam=size(Ind);

probability=zeros(tam2(2),1);
pval=probability;
mdia=probability;

xdata=zeros(numero,1);
ydata=xdata;
liminf=zeros(tam(2),1);
limsup=liminf;

for step=numero+1:1:tam(2)-numero-1
   for vez=1:numero
           ydata(vez,1)=double(ref(Ind(step-numero-1+vez)));
           y2data(vez,1)=double(ref(Ind(step+vez)));
           y3data(vez,1)=double(ref2(Ind(step+vez)));
           xdata(vez,1)=double(Ind(1,step-numero-1+vez))*1.000;
   end
       
   
    [p,S]=polyfit(xdata',ydata',grau);

    yfit = polyval(p,xdata);
    [Y,DELTA] = polyconf(p,xdata,S,'alpha',alpha, 'predopt', 'curve');
   
    liminf(step)=Y(numero)-DELTA(numero);
    limsup(step)=Y(numero)+DELTA(numero);
        if ref(Ind(step))<(Y(numero)-DELTA(numero))
            %alerta(Ind(step))=1;
             [h prob]=ttest2(ydata,y2data);
             
            if h==1
%                 conta=0;
%                 for nemaoir=1:numero
%                     if ref(Ind(step+nemaoir))<(Y(numero)-DELTA(numero))
%                     conta=conta+1;
%                     end
%                 end
%                 conta=conta/numero;              
%                 probability(Ind(step))=conta;
                pval(Ind(step))=1-prob;
                mdia(Ind(step))=mean(y3data);
            else
%                 probability(Ind(step))=0;
                pval(Ind(step))=0;
                mdia(Ind(step))=0;
            end            
        end
        
        
end
%[pval mdia]
%for step=numero+1:1:tam(2)-numero-1
pval(find(diajul<intdias(1)))=0;
pval(find(diajul>intdias(2)))=0;


if sum(pval)==0
    diamax=0;
    valor=0;
    valor_dist=0;    
else
    diamax=diajul(find(pval==max(pval)));
    valor=max(pval);
    valor_dist=mdia(find(pval==max(pval)));
end

if tam(2)>numero*2.5
    dj_min=diajul(Ind(numero+1));
    dj_max=diajul(Ind(tam(2)-numero-1));
end


end


