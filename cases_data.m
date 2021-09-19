function [Cases_60_plus,Cases_all]=cases_data(x)
%x is the name of the relevant mat file
load(x)


N_ill1=N_ill(:,1);
N_ill2=N_ill(:,2);

Who_Vec1=Who_Vec(:,1);
Who_Vec2=Who_Vec(:,2);

Who_ill1=Who_ill(:,1);
Who_ill2=Who_ill(:,2);

%load('Who_vec.mat')


    temp=0;
    temp1=0;
    temp2=0;
    clear Y1_all Y2_all Y3_all Y4_all Y5_all Y6_all
    clear Y1_ill Y2_ill Y3_ill Y4_ill Y5_ill Y6_ill
    clear Y1_ill_vac3 Y2_ill_vac3 Y3_ill_vac3 Y4_ill_vac3 Y5_ill_vac3 Y6_ill_vac3
    clear Y1_ill_vac2 Y2_ill_vac2 Y3_ill_vac2 Y4_ill_vac2 Y5_ill_vac2 Y6_ill_vac2
    clear Y1_ill_vac1 Y2_ill_vac1 Y3_ill_vac1 Y4_ill_vac1 Y5_ill_vac1 Y6_ill_vac1
    clear y1_ill y2_ill y_ill_young
    clear y1_ill_vac y2_ill_vac y_ill_vac_young y1_severe_vac y2_severe_vac y_severe_vac_young y1_severe_all y2_severe_all
    
    
    
    
    
    clear ill_vac
    clear ill
    clear All_ill
    clear ill_vac
    clear ill
    clear All_ill
    temp=0;
    temp1=0;
    temp2=0;
    [y,x]=hist(reshape(Who_ill2(:,:),1,size(Who_ill2,1)*size(Who_ill2,2)),0:299);
    RR=8395/max(smooth(y(2:70),0.1));
    %%
    
    for i=1:ll
        
        who_ill=reshape(Who_ill1(i,:),size(Who_ill1(i,:)));
        who_ill=who_ill(1:find(who_ill==0,1));
        
        All_ill(temp2+1:temp2+length(who_ill),1)=who_ill;
        All_ill(temp2+1:temp2+length(who_ill),2)=Who_ill2(i,1:length(who_ill));
        
        [val,pos]=intersect(who_ill,Who_Vec1(i,:));
        ill_vac(temp+1:temp+length(val),1)=val;% who is ill and vaccinted
        ill_vac(temp+1:temp+length(pos),2)=Who_ill2(i,pos);%when was infected
        ill_1=setdiff(who_ill,who_ill(pos));%ill_without vaccion
        [val_1,pos_1]=intersect(who_ill,ill_1);
        ill(temp1+1:temp1+length(val_1),1)=val_1;% who is ill and not vaccinted
        ill(temp1+1:temp1+length(pos_1),2)=Who_ill2(i,pos_1);%when was infected
        
        
        
        [val,pos]=intersect(Who_Vec1(i,:),who_ill);
        ill_vac(temp+1:temp+length(pos),3)=Who_Vec2(i,pos);
        temp=max(size(ill_vac,1)-1,0);
        temp1=max(size(ill,1)-1,0);
        temp2=max(size(All_ill,1)-1,0);
    end
    [Y1_all,Y2_all,Y3_all,Y4_all,Y5_all,Y6_all]=make_hist_ages_ill_2(All_ill);
    [Y1_ill,Y2_ill,Y3_ill,Y4_ill,Y5_ill,Y6_ill]=make_hist_ages_ill_2(ill);
    %%
    
    %Y1_all(1,:)=y1_all;Y2_all(1,:)=y2_all;Y3_all(1,:)=y3_all;Y4_all(1,:)=y4_all;Y5_all(1,:)=y5_all;Y6_all(1,:)=y6_all;
    %Y1_ill(1,:)=y1_ill;Y2_ill(1,:)=y2_ill;Y3_ill(1,:)=y3_ill;Y4_ill(1,:)=y4_ill;Y5_ill(1,:)=y5_ill;Y6_ill(1,:)=y6_ill;
    
    T=1:-0.05:0.05;
    for i=1:20
        for j=1:20
            
            [Y1_ill_vac,Y2_ill_vac,Y3_ill_vac,Y4_ill_vac,Y5_ill_vac,Y6_ill_vac]=make_hist_ages_vac_find2(ill_vac,1,T(i),T(j));
            
            
            
            
            clear y1_all y2_all
            y1_all(1,1:length(Y1_all))=smooth(Y1_all+Y2_all,0.1);
            y1_severe_all(1,:)=(y1_all(1,:)).*SM_fit(y1_all(1,:));
            
            
            y2_all(1,:)=smooth(Y3_all+Y4_all+Y5_all+Y6_all,0.1);
            y2_severe_all(1,:)=(y2_all(1,:)).*SM_fit1(y2_all(1,:));
            
            %%
            
            
            %%
            
            
            %%
            
            y1_ill(1,1:length(Y1_ill))=smooth(Y1_ill,0.1);
            y2_ill(1,1:length(Y1_ill))=smooth(Y2_ill+Y3_ill+Y4_ill+Y3_ill+Y5_ill+Y6_ill,0.1);
            y_young_ill(1,1:length(Y1_ill))=smooth(Y6_ill,0.1);
            
            %                     y_ill_over_50(1,1:length(Y1_ill))=smooth(Y3_ill,0.1)*RR;
            %                     y_ill_over_40(1,1:length(Y1_ill))=smooth(Y4_ill,0.1)*RR;
            %                     y_ill_under_40(1,1:length(Y1_ill))=smooth(Y5_ill,0.1)*RR;
            %                     y_ill_young(1,1:length(Y1_ill))=smooth(Y6_ill,0.1)*RR;
            %%
            y1_ill_vac(1,1:length(Y1_ill_vac))=smooth(Y1_ill_vac,0.1);
            y2_ill_vac(1,1:length(Y2_ill_vac))=smooth(Y2_ill_vac+Y3_ill_vac+Y4_ill_vac+Y5_ill_vac,0.1);
            %y_young_ill_vac(1,1:length(Y1_ill_vac))=smooth(Y6_ill_vac,0.1);
            length(y1_ill_vac)
            length(y1_ill)
            
            
            %                     y_ill_vac_over_50(1,:)=smooth(Y3_ill_vac,0.1)*RR;
            %                     y_ill_vac_over_40(1,:)=smooth(Y4_ill_vac,0.1)*RR;
            %                     y_ill_vac_under_40(1,:)=smooth(Y5_ill_vac,0.1)*RR;
            %                     y_ill_vac_young(1,:)=smooth(Y6_ill_vac,0.1)*RR;
            %                     %%
            temp=min(length(y1_ill_vac),length(y1_ill));
            Temp=(y1_ill(1:100)+y2_ill(1:100)+y1_ill_vac(1:100)+y2_ill_vac(1:100));
            RR=8348/max(Temp(1:100));
            y1_severe_vac(1,1:temp)=(RR*(y1_ill_vac(1,1:temp)+y1_ill(1,1:temp))).*SM_fit(RR*(y1_ill_vac(1,1:temp)+y1_ill(1,1:temp)));
            y2_severe_vac(1,1:temp)=(RR*(y2_ill_vac(1,1:temp)+y2_ill(1,1:temp))).*SM_fit1(RR*(y2_ill_vac(1,1:temp)+y2_ill(1,1:temp)));
            
            %                     y_severe_vac_50(1,:)=(y_ill_vac_over_50(1,:)+y_ill_over_50(1,:)).*fittedmodel_50(y_ill_vac_over_50(1,:)+y_ill_over_50(1,:))';
            %                     y_severe_vac_40(1,:)=(y_ill_vac_over_40(1,:)+y_ill_over_40(1,:)).*fittedmodel_40(y_ill_vac_over_40(1,:)+y_ill_over_40(1,:))';
            %                     y_severe_vac_under_40(1,:)=(y_ill_vac_under_40(1,:)+y_ill_under_40(1,:)).*fittedmodel_under_40(y_ill_vac_under_40(1,:)+y_ill_under_40(1,:))';
            %                     y_severe_vac_young(1,:)=(y_ill_vac_young(1,:)+y_ill_young(1,:)).*fittedmodel_under_40(y_ill_vac_young(1,:)+y_ill_young(1,:))';
            %
            %%
            Cases_60_plus(1:temp,i,j)=RR*(y1_ill_vac(1:temp)+y1_ill(1:temp));
            Cases_all(1:temp,i,j)=RR*(y2_ill(1:temp)+y2_ill_vac(1:temp)+y1_ill_vac(1:temp)+y1_ill(1:temp));
            
            %   end
            %end
        end
    end
end



