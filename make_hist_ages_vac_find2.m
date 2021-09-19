function[Y_60,Y_50,Y_40,Y_under_40,Y_young,Y_kids]= make_hist_ages_vac_find2(ill_vac,a1,a2,a3)

N1=5000:2000:15000;
clear new_case_40
clear new_case_50
clear new_case_60
clear new_case_70
clear new_case_young

clear New_case_40
clear New_case_50
clear New_case_60
clear New_case_70
clear New_case_young


ages=[100.0000   84.0145   74.9636   63.2161    36.4333   20.0726 0]/100;

ind=find(ill_vac(:,1)>=ages(2)*N1(4)&ill_vac(:,1)<ages(1)*N1(4));
y_60(1:length(ind),1)=ill_vac(ind,2);
y_60(1:length(ind),2)=being_ill_find(ill_vac(ind,2)-ill_vac(ind,3),a1,a2,a3);
y_60(1:length(ind),2)=floor(y_60(:,2)+rand(length(y_60(:,2)),1));


ind=find(ill_vac(:,1)>=ages(3)*N1(4)&ill_vac(:,1)<ages(2)*N1(4));
y_50(1:length(ind),1)=ill_vac(ind,2);
y_50(1:length(ind),2)=being_ill_find(ill_vac(ind,2)-ill_vac(ind,3),a1,a2,a3);
y_50(1:length(ind),2)=floor(y_50(:,2)+rand(length(y_50(:,2)),1));

ind=find(ill_vac(:,1)>=ages(4)*N1(4)&ill_vac(:,1)<ages(3)*N1(4));
y_40(1:length(ind),1)=ill_vac(ind,2);
y_40(1:length(ind),2)=being_ill_find(ill_vac(ind,2)-ill_vac(ind,3),a1,a2,a3);
y_40(1:length(ind),2)=floor(y_40(:,2)+rand(length(y_40(:,2)),1));

ind=find(ill_vac(:,1)>=ages(5)*N1(4)&ill_vac(:,1)<ages(4)*N1(4));
y_under_40(1:length(ind),1)=ill_vac(ind,2);
y_under_40(1:length(ind),2)=being_ill_find(ill_vac(ind,2)-ill_vac(ind,3),a1,a2,a3);
y_under_40(1:length(ind),2)=floor(y_under_40(:,2)+rand(length(y_under_40(:,2)),1));

ind=find(ill_vac(:,1)>=ages(6)*N1(4)&ill_vac(:,1)<ages(5)*N1(4));
y_young(1:length(ind),1)=ill_vac(ind,2);
y_young(1:length(ind),2)=being_ill_find(ill_vac(ind,2)-ill_vac(ind,3),a1,a2,a3);
y_young(1:length(ind),2)=floor(y_young(:,2)+rand(length(y_young(:,2)),1));



% 
Li=max(max(y_young(:,1)))-1;

[y,x]=hist(y_young(:,1).*y_young(:,2),0:Li);
Y_young=smooth(y(2:end),.1);

[y,x]=hist(y_60(:,1).*y_60(:,2),0:Li);
Y_60=smooth(y(2:end),.1);

[y,x]=hist(y_50(:,1).*y_50(:,2),0:Li);
Y_50=smooth(y(2:end),.1);

[y,x]=hist(y_40(:,1).*y_40(:,2),0:Li);
Y_40=smooth(y(2:end),.1);

[y,x]=hist(y_under_40(:,1).*y_under_40(:,2),0:Li);
Y_under_40=smooth(y(2:end),.1);

Y_kids=0;
%end
%