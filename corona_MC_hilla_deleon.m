function [N_ill,Who_Vec,Who_ill,Who_not_ill,Tt]=corona_MC_hilla_deleon(N,R,a1_trans,a1_infect,a2_infect)
ratio=(1.25./sqrt(R));
ages=[100,85,75,63,36];
load('Vac.mat')



N_60=round(N*ages(1)/100):-1:round(N*ages(2)/100);
N_60=N_60(randperm(length(N_60)));
%N_60=N_60(1:0.8*round(length(N_60)));
L_60=length(N_60);

N_50=round(N*ages(2)/100):-1:round(N*ages(3)/100);
N_50=N_50(randperm(length(N_50)));
%N_50=N_50(1:0.8*round(length(N_50)));
L_50=length(N_50);

N_40=round(N*ages(3)/100):-1:round(N*ages(4)/100);
N_40=N_40(randperm(length(N_40)));
%N_40=N_40(1:0.8*round(length(N_40)));
L_40=length(N_40);

N_young=round(N*ages(4)/100):-1:round(N*ages(5)/100);
N_young=N_young(randperm(length(N_young)));
%N_young=N_young(1:0.8*round(length(N_young)));
L_young=length(N_young);


D=100./(12:-1:1);
x0=linspace(0,5e3,5e3+1);
y0=linspace(0,5e3,5e3+1);
r0=sqrt(x0.^2+y0.^2);
patient_zero_x=x0(ceil(2e3*rand(1,1)));
patient_zero_y=y0(ceil(2e3*rand(1,1)));
days=[];
vec_days=0;
number_of_people=ceil(rand(1)*5);
n_not_ill=N-1;
n_ill=N-n_not_ill;
M_ill_x=patient_zero_x;
M_ill_y=patient_zero_y;

Who_free=(1:N);
Who_free=Who_free(randperm(length(Who_free)));
Who_free=Who_free(1:0.95*round(length(Who_free)));

Who_ill=[ceil(rand(1,1)*N),0];
Who_Vec=[];
Who_rec=[];
Who_not_ill=[Who_free,Who_rec];
Who_not_ill=Who_not_ill(Who_not_ill~=Who_ill(1,1));

n_not_ill=length(Who_not_ill);

M_not_ill_x(1:n_not_ill)=x0(ceil(2e3*rand(1,n_not_ill)));
M_not_ill_y(1:n_not_ill)=y0(ceil(2e3*rand(1,n_not_ill)));


%%
N_ill=0;

cycle = [zeros(1,21),0.9*ones(1,21),zeros(1,7),0.9*ones(1,201)]';
%j=6;
%old_people=Who_not_ill(Who_not_ill>groups(j-1) & Who_not_ill<=(groups(j)-N*0.15*0.1));
r_end=2e3;
%ratio=Ratio_09(:,3);
for n=1:300
    n;
    if n>25
       
        r_end=ratio(n)*2e3;
        if r_end>4e3
            r_end
            n
        end
         
    end
    
    
    clear temp
    clear temp2
    temp=ceil(normrnd(14,4,1,length(Who_ill(:,1)))/2+...
        normrnd(18,6,1,length(Who_ill(:,1)))/2);%time to recover
    %a1=
    a1=Who_ill((n-Who_ill(:,2)-temp')>0,1);
    a2=Who_ill((n-Who_ill(:,2)-temp')>0,2);
    
    Who_rec(end+1:end+length(a1),1)=a1;
    Who_rec(end+1-length(a2):end,2)=a2;
    active_case=setdiff(Who_ill(:,1),Who_rec(:,1));%numer of active cases
    N_ill(n,1)=length(setdiff(Who_ill(:,1),Who_rec(:,1)));
    N_ill(n,2)=length(find(active_case>N*0.85));
    
    if n>21&&n<98
        n;
        n_vaccine_60 =round(min(floor(L_60*(vac_60(n-20)-vac_60(n-21))),length(N_60)));
        n_vaccine_50 =round(min(floor(L_50*(vac_50(n-20)-vac_50(n-21))),length(N_50)));
        n_vaccine_40 =round(min(floor(L_40*(vac_40(n-20)-vac_40(n-21))),length(N_40)));
        n_vaccine_young =round(min(floor(L_young*(vac_young(n-20)-vac_young(n-21))),length(N_young)));
        
        Who_Vec(end+1:end+(n_vaccine_60),1)=N_60(1:n_vaccine_60);
        Who_Vec(end-n_vaccine_60+1:end,2)=n;
        N_60=N_60(n_vaccine_60+1:end);
        
        Who_Vec(end+1:end+(n_vaccine_50),1)=N_50(1:n_vaccine_50);
        Who_Vec(end-n_vaccine_50+1:end,2)=n;
        N_50=N_50(n_vaccine_50+1:end);
        
        Who_Vec(end+1:end+(n_vaccine_40),1)=N_40(1:n_vaccine_40);
        Who_Vec(end-n_vaccine_40+1:end,2)=n;
        N_40=N_40(n_vaccine_40+1:end);
        
        Who_Vec(end+1:end+(n_vaccine_young),1)=N_young(1:n_vaccine_young);
        Who_Vec(end-n_vaccine_young+1:end,2)=n;
        N_young=N_young(n_vaccine_young+1:end);
        Tt(n,1)=n_vaccine_60;
        Tt(n,2)=n_vaccine_50;
        Tt(n,3)=n_vaccine_40;
        Tt(n,4)=n_vaccine_young;
    end
    
    
    if n>97&&n<127
        n;
        %n_vaccine_60 =round(min(floor(L_60*(vac_60(end)-vac_60(end-1))),length(N_60)));
        n_vaccine_50 =round(min(floor(L_50*(vac_50(end)-vac_50(end-1))),length(N_50)));
        n_vaccine_40 =round(min(floor(L_40*(vac_40(end)-vac_40(end-1))),length(N_40)));
        n_vaccine_young =round(min(floor(L_young*(vac_young(end)-vac_young(end-1))),length(N_young)));
        
        Who_Vec(end+1:end+(n_vaccine_60),1)=N_60(1:n_vaccine_60);
        Who_Vec(end-n_vaccine_60+1:end,2)=n;
        N_60=N_60(n_vaccine_60+1:end);
        
        Who_Vec(end+1:end+(n_vaccine_50),1)=N_50(1:n_vaccine_50);
        Who_Vec(end-n_vaccine_50+1:end,2)=n;
        N_50=N_50(n_vaccine_50+1:end);
        
        Who_Vec(end+1:end+(n_vaccine_40),1)=N_40(1:n_vaccine_40);
        Who_Vec(end-n_vaccine_40+1:end,2)=n;
        N_40=N_40(n_vaccine_40+1:end);
        
        Who_Vec(end+1:end+(n_vaccine_young),1)=N_young(1:n_vaccine_young);
        Who_Vec(end-n_vaccine_young+1:end,2)=n;
        N_young=N_young(n_vaccine_young+1:end);
        Tt(n,1)=n_vaccine_60;
        Tt(n,2)=n_vaccine_50;
        Tt(n,3)=n_vaccine_40;
        Tt(n,4)=n_vaccine_young;
    end
    
    
    
    
    
    %M_free_x=M_free_x(1:end-n_vaccine);
    %M_free_y=M_free_y(1:end-n_vaccine);
    
    
    
    %
    %     if n>151 % reinfection after 150 days:
    %
    %         %save('M_free.mat','M_free_x','M_free_y','n','N_rec','x0','y0')
    %         X1=x0(ceil(2e3*rand(1,Who_rec(n-150,1))));
    %         Y1=y0(ceil(2e3*rand(1,N_rec(n-150,1))));
    %         for h=1:length (X1)
    %             LL=length(M_not_ill_x);
    %             M_not_ill_x(LL+1)=X1(h);
    %             M_free_y(LL+1)=Y1(h);
    %         end
    %
    %
    %
    %
    %     end
    %
    
    
    %n_ill=length(M_ill_x)
    
    
    l=3;
    %add new sick people - chance of 20%
    if ratio(n)<1.3&&n<100
        ch=0.8;
    else
        ch=0.95;
    end
    if rand(1,1)>ch
        n;
        M_ill_x(length(M_ill_x)+1)=x0(ceil(r_end*rand(1,1)));
        M_ill_y(length(M_ill_y)+1)=y0(ceil(r_end*rand(1,1)));
        
        M_not_ill_x=M_not_ill_x(1:end-1);
        M_not_ill_y=M_not_ill_y(1:end-1);
        
        Who_ill(end+1,1)=Who_not_ill(ceil(rand(1,1)*length(Who_not_ill)));
        Who_ill(end,2)=n;
        Who_not_ill=Who_not_ill(Who_not_ill~=Who_ill(end,1));
    end
    
    %M_free_x=M_free_x(1:end-n_vaccine);
    %M_free_y=M_free_y(1:end-n_vaccine);
    
    %%
    n_not_ill=size(Who_not_ill,2);
    x1_not_ill=normrnd(0,2*0.5*(10^l),1,n_not_ill);
    y1_not_ill=normrnd(0,2*0.5*(10^l),1,n_not_ill);
    
    n_ill=size(Who_ill,1);
    
    x1_ill=normrnd(0,2*0.5*(10^l),1,n_ill);
    y1_ill=normrnd(0,2*0.5*(10^l),1,n_ill);
    
    temp=normrnd(5,1,1,length(Who_ill(:,2)));
    temp(temp<0)=0;
    
    temp2=n-Who_ill(:,2);%number of days from infection
    %temp3=find(rand(1,length(days))>0.25);%people with symptoms
    temp4=find(temp<temp2);%sick people that are still sick.
    
    
    %                 x1_ill(intersect(temp3,temp4))=0;
    %                 y1_ill(intersect(temp3,temp4))=0;
    %
    M_ill_x=M_ill_x+x1_ill;
    M_ill_y=M_ill_y+y1_ill;
    M_ill_x=mod(M_ill_x,r_end);M_ill_y=mod(M_ill_y,r_end);
    
    
    %%
    
    
    M_not_ill_x=M_not_ill_x+x1_not_ill;
    M_not_ill_y=M_not_ill_y+y1_not_ill;
    
    
    M_not_ill_x=mod(M_not_ill_x,r_end);
    M_not_ill_y=mod(M_not_ill_y,r_end);
    %print(M_free_y(M_free_y>2e3))
    %print(M_free_y(M_free_x>2e3))
    
    
    t_M_ill=length(M_ill_x);
    if t_M_ill>0
        M_ill_x1=M_ill_x'*ones(1,n_not_ill);
        M_ill_y1=M_ill_y'*ones(1,n_not_ill);
        
        
        if n_not_ill>0
            M_not_ill_x1=ones(n_ill,1)*M_not_ill_x;
            M_not_ill_y1=ones(n_ill,1)*M_not_ill_y;
            
            T_x=abs(M_not_ill_x1-M_ill_x1);
            T_y=abs(M_not_ill_y1-M_ill_y1);
            T=sqrt(T_x.^2+T_y.^2);
            %temp5=intersect(temp3,temp4);%people with symptoms after ~5 days
            
%             map=ones(1,length(Who_ill));
%             
%             
%             if ~isempty(Who_Vec>0)
%                 
%                 map=ones(1,length(Who_ill));
%                 temp1=find(n-Who_Vec(:,2)>10);
%                 [val,pos]=intersect(Who_ill,Who_Vec(temp1,1));
%                 map(pos)=0.5;
%             end

            
            x_n=n-Who_ill(:,2);
            r_t=+exp(-(x_n-4).^2/2)+exp(-(x_n-5).^2/2)+exp(-(x_n-6).^2/2)+exp(-(x_n-7).^2/2);%chance of being infection a function of the day.
            r_t(r_t>1)=1;
            map=ones(size(r_t));
            size(map);
            if ~isempty(Who_Vec>0)
                
                temp1=find(n-Who_Vec(:,2)>10);
                [val,pos]=intersect(Who_ill,Who_Vec(temp1,1));
                size(pos);
                map(pos',:)=a1_trans;
            end
            n;
           if size(r_t,1)~=size(map,1)
               map=map';
           end
            size(r_t);
            size(map);
            
            numeber_of_oc=ceil(rand(1,t_M_ill)*2);
            a_symptom=round(rand(1,t_M_ill))/2+0.5;%half of the people are asymptomatic with P=1/2;
            %                         n;
            T(T<3)=3; %keepind SD
            temp3=find(a_symptom==1);
            temp5=intersect(temp3,temp4);%people with symptoms after ~5 days
            
            %                         for jj=1:length(temp5)
            %                             a=T(temp5(jj),:);
            %                             a(a<8)=8;
            %                             T(temp5(jj),:)=a;
            %                             min(a);
            %                         end
            t=1;
            P1=exp(-T.^2/2/2.4^2).*r_t.*map*number_of_people.*numeber_of_oc';
            
            map=ones(1,length(Who_not_ill));
            
            
            if ~isempty(Who_Vec>0)
               
                
                temp1=find((n-Who_Vec(:,2)>14));
                [val,pos]=intersect(Who_not_ill,Who_Vec(temp1,1));
                map(pos)=a2_infect-0.2;
                
                temp2=find((n-Who_Vec(:,2))>10 & (n-Who_Vec(:,2))<=14);
                [val,pos]=intersect(Who_not_ill,Who_Vec(temp2,1));
                map(pos)=a1_infect;
                
                temp3=find((n-Who_Vec(:,2))>7 & (n-Who_Vec(:,2))<=10);
                [val,pos]=intersect(Who_not_ill,Who_Vec(temp3,1));
                map(pos)=1;
            end
            P=sum(0.7*P1,1)+rand(1,length(Who_not_ill));
            P=floor(P.*map);
            P(P>1)=1;
            
            new_M_not_ill_x=M_not_ill_x(P<1);
            new_M_not_ill_y=M_not_ill_y(P<1);
            j;
            sum(P);
            if sum(P>0)
                M_ill_x(end+1:end+sum(P))=M_not_ill_x(P>0);
                M_ill_y(end+1:end+sum(P))=M_not_ill_y(P>0);
                Who_ill(end+1:end+sum(P),1)=Who_not_ill(P>0);
                Who_ill(end+1-sum(P):end,2)=n;
                Who_not_ill=setdiff(Who_not_ill,Who_ill(:,1));
            end
            n_not_ill=length(new_M_not_ill_x);
            M_not_ill_x=new_M_not_ill_x;
            M_not_ill_y=new_M_not_ill_y;
            %adding reinfection after 150 days:
            
            
        end
    end
    
    
end
%end

end