function [diamax valor valor_dist burn difdias twowayprob]=classifica_multiplo(sat, site, lcover, nir_ref, mask_gl, diajul);
% sat=1;
% site=1;
% nir_ref=medias_r3v;
% mask_gl=mask_glv;
% diajul=diajul_vgt;
% lcover=9;
ver=0;
timeb=[1 32 60 91 121 152 182 213 244 274 305 335 max(diajul)];

intdias=clima(site);

alpha=0.05;
grau=1;
if sat==1
    numero=10;    
elseif sat==2   
    numero=6;
end

ref=nir_ref';
mask=mask_gl';
validos=ref(find(mask~=1));
Ind=find(mask~=1);

diamax=zeros(12,1);
valor=zeros(12,1);
valor_dist=zeros(12,1);
burn=zeros(12,1);
difdias=zeros(12,1);   
twowayprob=zeros(12,1);

tam2=size(ref);
tam=size(Ind);

probability=zeros(tam2(2),1);
pval=probability;
mdia=probability;
probabildade=probability;
marca=probability;
marca2=probability;

xdata=zeros(numero,1);
ydata=xdata;

liminf=zeros(tam(2),1);
limsup=liminf;

for step=numero+1:1:tam(2)-numero-1
    for vez=1:numero
        ydata(vez,1)=double(ref(Ind(step-numero-1+vez)));
        y2data(vez,1)=double(ref(Ind(step+vez)));
        xdata(vez,1)=double(Ind(1,step-numero-1+vez))*1.000;
    end
    [h prob]=ttest2(ydata,y2data);
    pval(Ind(step))=1-prob;
    mdia(Ind(step))=mean(y2data);
    marca2(Ind(step))=2;
    [p,S]=polyfit(xdata',ydata',grau);
    %[lcover mdia(Ind(step))]
    probabildade(Ind(step))=max(obterprodist(lcover, mdia(Ind(step))/100, 1, site));
    
  
    yfit = polyval(p,xdata);
    [Y,DELTA] = polyconf(p,xdata,S,'alpha',alpha, 'predopt', 'curve');
    
    liminf(step)=Y(numero)-DELTA(numero);
    limsup(step)=Y(numero)+DELTA(numero);
    
    if ref(Ind(step))<(Y(numero)-DELTA(numero))
        [h prob]=ttest2(ydata,y2data);
        
        if h==1
            marca(Ind(step))=1;
        else
            marca(Ind(step))=0;
        end
    end
    
    
end

probprod=pval.*probabildade;

agua=zeros(12,1);
aguatag=ones(12,1)*3;

if sum(probprod)>0.01 | lcover>0
    
%[probprod marca marca2]

superiores=numel(find(marca==1));
if superiores>0
arde=find(probprod==max(probprod(marca==1)));
if numel(arde)==1
    marca2(arde)=1;
end
end



marca2(find(diajul<intdias(1)))=2;
marca2(find(diajul>intdias(2)))=2;
% 
% burn=diajul(find(pval==max(pval) & marca==1));

% serie dos p-val significativos e do valor da reflectancia pos 
% [pval e mdia]

mes=0;
for i=2:13
    mes=mes+1;
    I=find(diajul>=timeb(i-1) &  diajul<timeb(i));
    %numel(I)
    if numel(I)>0
        
        %tag=find(pval==max(pval(I)))
        tag=find(probprod==max(probprod(I)));

        valor_dist(mes)=mdia(tag(1));
        valor(mes)=pval(tag(1));
        diamax(mes)=diajul(tag(1));
        burn(mes)=marca2(tag(1)); 
        if(tag(1)==1)
            difdias(mes)=0;
        else
            difdias(mes)=diajul(tag(1))-diajul(tag(1)-1);
        end
        twowayprob(mes)=probprod(tag(1));
    else
        valor_dist(mes)=0;
        valor(mes)=0;
        diamax(mes)=0;
        burn(mes)=0;
        difdias(mes)=0;
        twowayprob(mes)=0;
   end
end

else
    
    valor_dist=agua;
    valor=agua;
    diamax=agua;
    burn=aguatag;
    difdias=agua;
    twowayprob=agua;
    
    
end

if ver==1
figure; 


subplot1 = subplot(5,1,1,...
    'XTickLabel',{'Jan' 'Fev' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Ago' 'Sep' 'Oct' ' Nov' 'Dec'},...
    'XTick',[1 32 60 91 121 152 182 213 244 274 305 335]);
box(subplot1,'on');
hold(subplot1,'all');


plot(diajul_vgt,medias_r3v)
hold on; plot(diajul_vgt(find(mask_glv==0)), medias_r3v(find(mask_glv==0)), 'r')
xlim([0 367])
ylabel('ref')


subplot2 = subplot(5,1,2,...
    'XTickLabel',{'Jan' 'Fev' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Ago' 'Sep' 'Oct' ' Nov' 'Dec'},...
    'XTick',[1 32 60 91 121 152 182 213 244 274 305 335]);
box(subplot2,'on');
hold(subplot2,'all');


plot(diajul_vgt, pval); hold on; plot(diajul_vgt(find(mask_glv==0)), pval(find(mask_glv==0)), 'r')
hold on; plot(diamax, valor, 'ok')
hold on; plot(diajul_vgt(find(mask_glv==0)), probabildade(find(mask_glv==0)), 'k')
hold on; plot(diajul_vgt(find(mask_glv==0)), probabildade(find(mask_glv==0)).*pval(find(mask_glv==0)), 'g')
xlim([0 367])
ylabel('p-val')



subplot4 = subplot(5,1,3,...
    'XTickLabel',{'Jan' 'Fev' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Ago' 'Sep' 'Oct' ' Nov' 'Dec'},...
    'XTick',[1 32 60 91 121 152 182 213 244 274 305 335]);
box(subplot4,'on');
hold(subplot4,'all');

plot(diajul_vgt, mdia); hold on; plot(diajul_vgt(find(mask_glv==0)), mdia(find(mask_glv==0)), 'r')
hold on; plot(diamax, valor_dist, 'ok')
xlim([0 367])
ylabel('med.ref.pos')

subplot5 = subplot(5,1,4,...
    'XTickLabel',{'Jan' 'Fev' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Ago' 'Sep' 'Oct' ' Nov' 'Dec'},...
    'XTick',[1 32 60 91 121 152 182 213 244 274 305 335]);
box(subplot5,'on');
hold(subplot5,'all');

bar(diamax, burn, 'r')
xlim([0 367])
ylabel('med.ref.pos')

subplot6 = subplot(5,1,5,...
    'XTickLabel',{'Jan' 'Fev' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Ago' 'Sep' 'Oct' ' Nov' 'Dec'},...
    'XTick',[1 32 60 91 121 152 182 213 244 274 305 335]);
box(subplot6,'on');
hold(subplot6,'all');

bar(diamax, difdias, 'r')
xlim([0 367])
ylabel('med.ref.pos')
end


end

