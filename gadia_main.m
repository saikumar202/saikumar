

clc;
clear all
close all
Nu=2;
Nf=2;
pf=[10 1];
C=3;
B=1;
K=[3 5];
T=120;
lambda=[1,2; 3,4;];
Pt = 0:1:5;
pc = 0.1:0.1:0.9;

x=zeros(Nu,Nf);
count=zeros(1,Nu);
add_item=zeros(Nu,Nu);
index=ones(Nu,Nf);
add_max=zeros(Nf,2);  % maximum priority value of each file
add_sum=zeros(Nu,Nf); % priority value
offloading_ratio=0;
jx=ones(max(K)+1,1);
for (i=3:max(K)+1) jx(i)=jx(i-1)*(i-1);
end

% initialize the priority values
for (j=1:Nu)
    for (i=1:Nu)
        add_item(j,i)=1/Nu*(1-exp(-T*lambda(i,j)));
    end
end  
for (k=1:Nf)
    add_sum(:,k)=pf (k)/K(k)*add_item(:,:)*index(:,k);
    [add_max(k,1),add_max(k,2)]=max(add_sum(:,k));
end

% greedy algorithm
for (j=1:Nu*C)
    
    % select the element with maximum priority value
    [~,Nav]=max(add_max(:,1));
    Uav=add_max(Nav,2);
    offloading_ratio=offloading_ratio+add_max(Nav,1);
    
    % update the remaining set, which includes the remaining elements that can be added considering the constraint
    count(Uav)=count(Uav)+1;
    if (count(Uav)>=C)
        add_sum(Uav,:)=ones(1,Nf)*(-inf);
        for (k=1:Nf)
            [add_max(k,1),add_max(k,2)]=max(add_sum(:,k));
        end
    end    
    x(Uav,Nav)=x(Uav,Nav)+1;
    if (x(Uav,Nav)>=K(Nav))
        add_sum(Uav,Nav)=-inf;
    end
    
    % update the priority values
	for (j1=1:Nu)
		if (add_sum(j1,Nav)>=0) add_sum(j1,Nav)=0;
        end
    end  
    for (i=1:Nu)
		if (K(Nav)-x(i,Nav)-1>=0)
			temp_prob=[1 0.1 0.01 0.001];
            if (add_sum(i,Nav)>=0)
                add_sum(i,Nav)=add_sum(i,Nav)+1;
            end
            for (j1=1:Nu)
                if (i ~= j1 && add_sum(j1,Nav)>=0)
                    if (x(j1,Nav)==0) add_sum(j1,Nav)=add_sum(j1,Nav)+gammainc(T*lambda(i,j1),ceil((x(j1,Nav)+1)/B))*temp_prob(K(Nav)-x(i,Nav))*pf(Nav)/Nu/K(Nav);
                    else if (K(Nav)-x(i,Nav)-x(j1,Nav)-1>=0)
                            temp_k=K(Nav)-x(i,Nav)-x(j1,Nav);
                            temp_a=zeros(temp_k,temp_k+1);
                            prob_j0=exp(-T*lambda(i,j1));
                            for (li=1:temp_k) 
                                temp_a(li,li)=1;
                                temp_a(li,temp_k+1)=temp_prob(temp_k-li+1)/prob_j0;
                            end

                            prob_j=zeros(temp_k-1,1);
                            for (z=1:min(x(j1,Nav)-1,temp_k-1))
                                if (mod(z,B) == 0)
                                    prob_j(z)=(T*lambda(i,j1))^(z/B)/(jx(z/B+1));
                                else
                                    prob_j(z)=0;
                                end
                            end
                            if (x(j1,Nav)<=temp_k-1) prob_j(x(j1,Nav))=(1-gammainc(floor(x(j1,Nav)/B),T*lambda(i,j1)))/prob_j0;
                            end
                            for (li=1:temp_k-1)
                                for (lj=li+1:temp_k)
                                    temp_a(li,lj)=prob_j(lj-li);
                                end
                            end
                            for (li=2:temp_k) temp_a(1,:)=temp_a(1,:)-temp_a(1,li)*temp_a(li,:);
                            end
                            add_sum(j1,Nav)=add_sum(j1,Nav)+gammainc(T*lambda(i,j1),ceil((x(j1,Nav)+1)/B))*temp_a(1,temp_k+1)*pf(Nav)/Nu/K(Nav);
                        end
                    end
                end
            end
        end
    end
    [add_max(Nav,1),add_max(Nav,2)]=max(add_sum(:,Nav));

 [St,st]=throughput(add_max(Nav,1),add_max(Nav,2));
  
end



figure;
plot (Pt,St(:,1),'k-*','LineWidTh',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (Pt,St(:,2),'r-*','LineWidth',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (Pt,St(:,3),'g-*','LineWidTh',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (Pt,St(:,4),'k-o','LineWidth',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (Pt,St(:,5),'r-o','LineWidth',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (Pt,St(:,6),'g-o','LineWidth',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (Pt,St(:,13),'m-o','LineWidth',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (Pt,St(:,14),'m-*','LineWidth',2,'MarkerFaceColor','b','Markersize',7);
hold off
grid on
xlabel('transmot power of the ps');
ylabel('sum throughput bps/hz');
legend('proposed algorithum d=10m','tdma d=10m','oet d=10m','proposed algorithum d=8m','tdma d=8m','oet d=8m','GADIA d=10m','GADIA d=8m');
title(' throughput performance comparison of different schemes ')

figure;
plot (pc,st(:,7),'k-*','LineWidTh',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (pc,st(:,8),'r-*','LineWidth',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (pc,st(:,9),'g-*','LineWidTh',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (pc,st(:,10),'k-o','LineWidth',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (pc,st(:,11),'r-o','LineWidth',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (pc,st(:,12),'g-o','LineWidth',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (pc,st(:,15),'m-o','LineWidth',2,'MarkerFaceColor','b','Markersize',7);
hold on
plot (pc,st(:,16),'m-*','LineWidth',2,'MarkerFaceColor','b','Markersize',7);
hold off
grid on

xlabel(' circuit power consumption');
ylabel('sum throughput(bps/hz)');
legend('proposed algorithum N=6','tdma N=6','oet d=10m','proposed algorithum d=8m','tdma d=8m','oet d=8m','GADIA d=10m','GADIA d=8m');
title(' throughput performance comparison of different schemes ')





