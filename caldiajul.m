function diajul=caldiajul(mesfil,diafil)

tam=size(diafil);
diajul=zeros(tam(1),1);

for i=1:tam(1)

if mesfil(i)==1    
    diajul(i)=diafil(i);    
elseif mesfil(i)==2    
    diajul(i)=diafil(i)+31;
elseif mesfil(i)==3    
    diajul(i)=diafil(i)+31+28;
elseif mesfil(i)==4    
    diajul(i)=diafil(i)+31+28+31;
elseif mesfil(i)==5    
    diajul(i)=diafil(i)+31+28+31+30;
elseif mesfil(i)==6    
    diajul(i)=diafil(i)+31+28+31+30+31;
elseif mesfil(i)==7    
    diajul(i)=diafil(i)+31+28+31+30+31+30;
elseif mesfil(i)==8    
    diajul(i)=diafil(i)+31+28+31+30+31+30+31;
elseif mesfil(i)==9    
    diajul(i)=diafil(i)+31+28+31+30+31+30+31+31;
elseif mesfil(i)==10    
    diajul(i)=diafil(i)+31+28+31+30+31+30+31+31+30;
elseif mesfil(i)==11    
    diajul(i)=diafil(i)+31+28+31+30+31+30+31+31+30+31;
elseif mesfil(i)==12    
    diajul(i)=diafil(i)+31+28+31+30+31+30+31+31+30+31+30;
end

end

end