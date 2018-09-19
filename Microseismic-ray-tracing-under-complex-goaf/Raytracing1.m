function [rays]=Raytracing1(data,Received)
%%%%¹ì¼£ ±àºÅ
% received=511;
for jj=1:size(Received,2)
lu=1;received=Received(jj);XX=[];YY=[];ZZ=[];
while ~isnan(data(received).before)
    XX(lu)=data(received).x;YY(lu)=data(received).y;ZZ(lu)=data(received).z;lu=lu+1;
    received=data(received).before;
end
XX(lu)=data(received).x;YY(lu)=data(received).y;ZZ(lu)=data(received).z;
rays{jj}(:,1)=XX;rays{jj}(:,2)=YY;rays{jj}(:,3)=ZZ;
end
end