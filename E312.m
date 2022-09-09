% this code produces Fig. 2 and corresponding data in Section 3.1.2

clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following code calculates the regret of policy pi_1 in Example 2 
% (Section 3.1.2) for each parameter set (a,b) in [8,10]*[4,5]

acc=0.001; % grid length of x-axis
accp=0.0001; % search step of determining the optimal price
B=zeros(1+1.25/acc,1+2/acc); % regret of pi_1 for each parameter set (a,b)

for i=1:1.25/acc+1 % this loop finds the regret of pi_1 for the left region of [8,10]*[4,5]
    cx(i)=8;
    cy(i)=(i-1)/1.25*acc+4;
    clear mreg
    for j=1:0.6/accp+1 % for each possible price
        p=0.7+(j-1)*accp;
        clear reg
        for k=1:1.25/acc+2-i
            ex=cx(i)+(k-1)*acc;
            ey=cy(i)+(k-1)/1.25*acc;
            rev=ex^2/4/ey;
            reg(k)=rev-p*(ex-ey*p);
        end
        mreg(j)=max(reg); % find the largest regret along this line segment for this price
    end
    mtr=mreg(1);
    mtt=1;
    for j=1:0.6/accp+1
        if mtr>mreg(j)
            mtr=mreg(j);
            mtt=j;
        end
    end
    popt=0.7+(mtt-1)*accp; % the optimal price along this line segment
    for k=1:1.25/acc+2-i
        ex=cx(i)+(k-1)*acc;
        ey=cy(i)+(k-1)/1.25*acc;
        B(k+i-1,k)=ex^2/4/ey-popt*(ex-ey*popt);
    end     
end

for i=1:1.25/acc+1 % this loop finds the regret of pi_1 for the middle region of [8,10]*[4,5]
    cx(i)=10;
    cy(i)=5-(i-1)/1.25*acc;
    clear mreg
    for j=1:0.6/accp+1 % for each possible price
        p=0.7+(j-1)*accp;
        clear reg
        for k=1:1.25/acc+2-i
            ex=cx(i)-(k-1)*acc;
            ey=cy(i)-(k-1)/1.25*acc;
            rev=ex^2/4/ey;
            reg(k)=rev-p*(ex-ey*p);
        end
        mreg(j)=max(reg); % find the largest regret along this line segment for this price
    end
    mtr=mreg(1);
    mtt=1;
    for j=1:0.6/accp+1
        if mtr>mreg(j)
            mtr=mreg(j);
            mtt=j;
        end
    end
    popt=0.7+(mtt-1)*accp; % the optimal price along this line segment
    for k=1:1.25/acc+2-i
        ex=cx(i)-(k-1)*acc;
        ey=cy(i)-(k-1)/1.25*acc;        
        B(1.25/acc-i-k+3,2/acc+2-k)=ex^2/4/ey-popt*(ex-ey*popt);
    end        
end

for i=2:(2-1.25)/acc % this loop finds the regret of pi_1 for the right region of [8,10]*[4,5]
    cx(i)=8+(i-1)*acc;
    cy(i)=4;
    clear mreg
    for j=1:0.6/accp+1 % for each possible price
        p=0.7+(j-1)*accp;
        clear reg
        for k=1:1.25/acc+1
            ex=cx(i)+(k-1)*acc;
            ey=cy(i)+(k-1)/1.25*acc;
            rev=ex^2/4/ey;
            reg(k)=rev-p*(ex-ey*p);
        end
        mreg(j)=max(reg); % find the largest regret along this line segment for this price
    end
    mtr=mreg(1);
    mtt=1;
    for j=1:0.6/accp+1
        if mtr>mreg(j)
            mtr=mreg(j);
            mtt=j;
        end
    end
    popt=0.7+(mtt-1)*accp; % the optimal price along this line segment
    for k=1:1.25/acc+1
        ex=cx(i)+(k-1)*acc;
        ey=cy(i)+(k-1)/1.25*acc;        
        B(k,i+k-1)=ex^2/4/ey-popt*(ex-ey*popt);
    end        
end

% since the size of B is 1251*2001, we transform it 
% into 801*2001 by interpolation, which is the size of A

C=zeros(801,2001); % regret of pi_1 for each parameter set (a,b)
C(801,:)=B(1251,:);
for i=1:16:785
    j=(i-1)*1.25*1.25+1;
    C(i,:)=B(j,:);
    C(i+1,:)=7*B(j+1,:)/16+9*B(j+2,:)/16;
    C(i+2,:)=14*B(j+3,:)/16+2*B(j+4,:)/16;
    C(i+3,:)=5*B(j+4,:)/16+11*B(j+5,:)/16;
    C(i+4,:)=12*B(j+6,:)/16+4*B(j+7,:)/16;
    C(i+5,:)=3*B(j+7,:)/16+13*B(j+8,:)/16;
    C(i+6,:)=10*B(j+9,:)/16+6*B(j+10,:)/16;
    C(i+7,:)=1*B(j+10,:)/16+15*B(j+11,:)/16;
    C(i+8,:)=8*B(j+12,:)/16+8*B(j+13,:)/16;
    C(i+9,:)=15*B(j+14,:)/16+1*B(j+15,:)/16;
    C(i+10,:)=6*B(j+15,:)/16+10*B(j+16,:)/16;
    C(i+11,:)=13*B(j+17,:)/16+3*B(j+18,:)/16;
    C(i+12,:)=4*B(j+18,:)/16+12*B(j+19,:)/16;
    C(i+13,:)=11*B(j+20,:)/16+5*B(j+21,:)/16;
    C(i+14,:)=2*B(j+21,:)/16+14*B(j+22,:)/16;
    C(i+15,:)=9*B(j+23,:)/16+7*B(j+24,:)/16;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following code calculates the regret of policy pi_2 in Example 2 
% (Section 3.1.2) for each parameter set (a,b) in [8,10]*[4,5]

G=zeros(801,2001);
for i=1:801
    for j=1:2001
        ex=(j-1)/1000+8;
        ey=(i-1)/800+4;
        if ex-ey<4.9
            G(i,j)=ex^2/4/ey-1.09375*(ex-ey*1.09375);
        else
            G(i,j)=ex^2/4/ey-1*(ex-ey*1);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following code calculates the regret of policy pi_3 in Example 2 
% (Section 3.1.2) for each parameter set (a,b) in [8,10]*[4,5]

H=0.01*G+0.99*C;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following code calculates the regret of different policies in Example  
% 2 (Section 3.1.2) for each parameter set (a,b) in [8,10]*[4,5]
save E312C.dat C -ascii; % regret of policy pi_1
save E312G.dat G -ascii; % regret of policy pi_2
save E312H.dat H -ascii; % regret of policy pi_3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following code draws Fig. 2
close all

fig2=figure(2);
fig2.Position=[50 50 1350 450];
subplot(1,2,1)
hold on
contourf(C,'ShowText','on')
xticks([1 401 801 1201 1601 2001])
xticklabels({'8','8.4','8.8','9.2','9.6','10'})
yticks([1 161 321 481 641 801])
yticklabels({'4','4.2','4.4','4.6','4.8','5'})
plot(2001,801,'ro','MarkerSize',20)
plot(751,1,'ro','MarkerSize',20)
title('Total regret for \pi_1')
set(gca,'FontSize',26);
hold off

subplot(1,2,2)
contourf(G,'ShowText','on')
xticks([1 401 801 1201 1601 2001])
xticklabels({'8','8.4','8.8','9.2','9.6','10'})
yticks([1 161 321 481 641 801])
yticklabels({'4','4.2','4.4','4.6','4.8','5'})
title('Total regret for \pi_2')
set(gca,'FontSize',26);

saveas(fig2,'Figure2.png')