%%%%%点云数据的射线追踪法
%%输入参数  判断两点的连线是否穿过空区，多个空区的话，需要区分开，KD1...
%%%%单个空区情况
clear;clc;
data=[];V=4000;
%% 空区1  mian DD Y PD
mian=load('C:\Users\Administrator\Desktop\A new method of microseismic ray tracing under complex goaf model\单空区模型mian.txt');
DD=load('C:\Users\Administrator\Desktop\A new method of microseismic ray tracing under complex goaf model\单空区模型dd.txt');
% Y=load('C:\Users\Administrator\Desktop\A new method of microseismic ray tracing under complex goaf model\y.txt');

Y=funpd(mian,DD);
% point=input('输入体内部一点坐标:'); %%%%[4000,2100,-765]
point=[60,390,0];
PD=Y*[point';1];
XI=find(PD<0);EY=find(PD==0);DY=find(PD>0);%%%小于0为0，等于0为1，大于0为2
PD(XI)=0;PD(EY)=1;PD(DY)=2;
mian=mian+1;


%%
Xmax=max(DD(:,1));Xmin=min(DD(:,1));Ymax=max(DD(:,2));Ymin=min(DD(:,2));Zmax=max(DD(:,3));Zmin=min(DD(:,3));
%%%data中加入震源点与接收点
%%单空区模型 研究区域 80 40 30 空区大小40 40 30 
Source=[Xmin-20,Ymin,Zmin+20];
Receiver1=[Xmax+20,Ymin,Zmin];
for i=1:5
    for j=1:4
        Receiver((i-1)*4+j,1)=Receiver1(1); Receiver((i-1)*4+j,2)=Receiver1(2)+(i-1)*10; Receiver((i-1)*4+j,3)=Receiver1(3)+(j-1)*10;
    end
end
c=1;
for i=1:size(Source,1)
    data(c).x=Source(i,1);data(c).y=Source(i,2);data(c).z=Source(i,3);data(c).t=inf;data(c).before=NaN;data(c).flag1=0;data(c).flag2=0;c=c+1;
end
for i=1:size(Receiver,1)
    data(c).x=Receiver(i,1);data(c).y=Receiver(i,2);data(c).z=Receiver(i,3);data(c).t=inf;data(c).before=NaN;data(c).flag1=0;data(c).flag2=0;c=c+1;
end
received=2:size(Receiver,1)+1;
Snum=1;Jmin=1;
data(Snum).t=0;data(Snum).flag1=2;data(Snum).flag2=2;

%% 
[data]=Travetime1(data,DD,V,Snum,Y,PD,Jmin);

[rays]=Raytracing1(data,received);

%% 到时矩阵已经得到 可视化验证  路径轨迹
 figure(1);
for i=1:size(rays,2)
    X(i)=plot3(rays{i}(:,1),rays{i}(:,2),rays{i}(:,3),'r');hold on;
end
 for i=1:size(mian,1)
     for j=1:size(mian,2)-1
         Xmian(1,j)=DD(mian(i,j+1),1);Ymian(1,j)=DD(mian(i,j+1),2);Zmian(1,j)=DD(mian(i,j+1),3);
     end
 miand(i)=patch(Xmian,Ymian,Zmian,[0.5,0.3,0.5]);hold on;
 end
 dian=scatter3(DD(:,1),DD(:,2),DD(:,3),'og');

function Y=funpd(mian,DD)
mian=mian+1;
for i=1:size(mian,1)
point1=DD(mian(i,2),:);point2=DD(mian(i,3),:);point3=DD(mian(i,4),:);
syms x y z;
D=[ones(4,1),[[x,y,z];point1;point2;point3]];
 detd = det(D);str = char(detd);
    amark = find(str=='x');
    bmark = find(str=='y');
    cmark = find(str=='z');
    dd = double(coeffs(detd));
    %------------------------------------------- 
    if(size(amark,2))
        if size(dd,2)==1
            a = coeffs(detd,x);   a= double(a(1));
        else
            a = coeffs(detd,x);   a = double(a(2));%系数按照升幂顺序排列
        end
    else
        a = 0;
    end
%------------------------------------------- 
    if(size(bmark,2))
        if size(dd,2)==1
            b = coeffs(detd,y);   b= double(b(1));
        else
            b = coeffs(detd,y);   b = double(b(2));%系数按照升幂顺序排列
        end
    else
        b = 0;
    end
%------------------------------------------- 
    if(size(cmark,2))
        if size(dd,2)==1
            c = coeffs(detd,z);   c= double(c(1));
        else
            c = coeffs(detd,z);   c = double(c(2));%系数按照升幂顺序排列
        end
    else
        c = 0;
    end
    
%------------------------------------------- 
if size(dd,2)==4
    d=dd(1);
end
if size(dd,2)==3
    if a~=0&&b~=0&&c~=0
        d=0;
    else
        d=dd(1);
    end
end
if size(dd,2)==2
    if (a~=0&&b~=0)||(a~=0&&c~=0)||(c~=0&&b~=0)
        d=0;
    else
        d=dd(1);
    end
end
if size(dd,2)==1
    d=0;
end
%-------------------------------------------
Y(i,:)=[a,b,c,d];
end
end
