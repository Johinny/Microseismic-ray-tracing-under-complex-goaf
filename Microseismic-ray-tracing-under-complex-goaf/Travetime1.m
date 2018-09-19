function [data]=Travetime1(data,DD,V,Snum,Y,PD,Jmin)
%%%DD空区边界点（n×3)，mian空区表面模型面片(N×(n+1) 第一列为面的标记，
%%%%V速度，data出空区点外的其他离散点（为一结构体1×N，单个结构体属性x,y,z坐标值,到时t（未知inf，出震源外）,
%%%%before前一个点,flag1其值为0,1,2(普通点、边界点、震源点),flag2其值0,1,2(未知到时点 不确定到时点 确定到时点 )
%%%%输出data，得到属性t
% DD=DD0{1};DD1=DD0{2};Y=Y0{1};Y1=Y0{2};PD=PD0{1};PD1=PD0{2};
ccc=length(data);
for i=1:size(DD,1)
    data(ccc+i).x=DD(i,1);data(ccc+i).y=DD(i,2);data(ccc+i).z=DD(i,3);
    data(ccc+i).t=inf;data(ccc+i).before=NaN;data(ccc+i).flag1=1;data(ccc+i).flag2=0;
    KD(i)=ccc+i;
end

%% 计算开始
dd=1;ww=1;Wt=[];

for i=1:size(data,2)
    
    if i~=Snum
     judge0=Judgeby0([data(Snum).x,data(Snum).y,data(Snum).z],[data(i).x,data(i).y,data(i).z],Y,PD,Jmin);
    if judge0==0
        data(i).t=norm([data(Snum).x-data(i).x,data(Snum).y-data(i).y,data(Snum).z-data(i).z],2)/V;
        data(i).before=Snum;data(i).flag2=2;Dt(dd)=i;dd=dd+1;
    else
        Wt(ww)=i;ww=ww+1;
    end
    end
end
tic
 
 %% 开始循环  循环停止条件  未知点为Wt 0，不确定点为Ut 0
%  Ut=[];sc=1;&&isempty(Ut)

 while ~(isempty(Wt))
      cseise=intersect(KD,Dt); %%%已知到时且在边界上的点
         %%%%%生成已知到时边界点矩阵ccseise
        ddd=1; Dtt=[];
     for i=1:size(cseise,2)
         
         for ii=1:size(Wt,2) %%%%%给未知到时、与不确定到时点计算到时
             
%              judge0=Judgeby0([data(cseise(i)).x,data(cseise(i)).y,data(cseise(i)).z],[data(Wt(ii)).x,data(Wt(ii)).y,...
%                  data(Wt(ii)).z],mian,DD);
               judge0=Judgeby0([data(cseise(i)).x,data(cseise(i)).y,data(cseise(i)).z],[data(Wt(ii)).x,data(Wt(ii)).y,...
                 data(Wt(ii)).z],Y,PD,Jmin);
%                judge1=Judgeby0([data(cseise(i)).x,data(cseise(i)).y,data(cseise(i)).z],[data(Wt(ii)).x,data(Wt(ii)).y,...
%                  data(Wt(ii)).z],Y1,PD1,Jmin);
              if judge0==0
                 tt=data(cseise(i)).t+norm([data(cseise(i)).x-data(Wt(ii)).x,data(cseise(i)).y-data(Wt(ii)).y,...
                     data(cseise(i)).z-data(Wt(ii)).z],2)/V;
                 if tt<=data(Wt(ii)).t  %%%%
                     data(Wt(ii)).t=tt;
                     data(Wt(ii)).before=cseise(i);
                 end
                 
   
                 %%%%%判断是否为到时确定点
                 judge1=Judge1([data(Wt(ii)).x,data(Wt(ii)).y,data(Wt(ii)).z],[data(cseise(i)).x,data(cseise(i)).y,...
                     data(cseise(i)).z],cseise,data,V,cseise(i));
                 if judge1==1 %%%%确定到时的点   在未知点集合Wt中去掉
                     data(Wt(ii)).flag2=2; Dtt(ddd)=Wt(ii);ddd=ddd+1; %%%%记录到时确定的点 
                 else   %%%%不确定到时的点
%                      data(Wt(ii)).flag2=1;
%                      Ut(sc)=Wt(ii);sc=sc+1;
                 end
              else   %%%%未知到时点里面还不知道的点
                 
              end
              
         end
         if i==63
             
         end
     end
     Dtt=unique(Dtt);
     %%%%%%更新Wt
     if isempty(Dtt)
     else
     for iii=1:size(Dtt,2)
        Wt(find(Wt==Dtt(iii)))=[];
     end
     end
     %%%%%%更新Dt
     Dt=Dtt;
    
 end
 toc
%  data(1).before=NaN;
 
end

function judge1=Judge1(A,B,cseise,data,V,index) %%%判断A此点是否为确定到时点 通过到边界点的距离 B本边界点，C其他边界点
Db=norm([A(1)-B(1),A(2)-B(2),A(3)-B(3)],2)/V+data(index).t;
for i=1:size(cseise,2)
     Dc(i)=norm([A(1)-data(cseise(i)).x,A(2)-data(cseise(i)).y,A(3)-data(cseise(i)).z],2)/V+data(cseise(i)).t;
end
if Db<=min(Dc)
    judge1=1;
else
    judge1=0;
end
end


function judge0=Judgeby0(A,B,Y,PD,Jmin)
%%%%AB线段上生成点
judge0=0;
D=norm(A-B,2);c=1;
  
while D>Jmin
point=[];
 for i=1:2^(c-1)
     xx=D/(2^c);bil=(2*i-1)*xx/D;
     point(i,:)=bil*B+(1-bil)*A;

    PD1=funpd(Y,point(i,:));
     if PD==PD1
         
       judge0=1;break;
     end
 end
 
 c=c+1;D=D/2;
end

if D<=Jmin
     point=[];
     point(1,:)=(A+B)/2;
     PD1=funpd(Y,point(1,:));
    if PD==PD1
       judge0=1;
    end
end

end
function PD1=funpd(Y,point)
PD1=Y*[point';1];
E=10E-6;
XI=find(PD1<-E);EY=[find(PD1>=-E);find(PD1<=E)];DY=find(PD1>E);%%%小于0为0，等于0为1，大于0为2
PD1(XI)=0;PD1(EY)=1;PD1(DY)=2;
end