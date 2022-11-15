%% useless
fid = fopen('sullivan.xyz','r');
fread(fid,5,'integer*4');
num = fread(fid,1,'integer*4');
xyz = fread(fid,num/8,'real*8');
fclose(fid);

%%
% a=37/180*pi;
% b=128/180*pi;
% c=173/180*pi;
tt=0.5;
a=0.35;
b=2.25;
c=1.5;
syms t
qa=[1 0 0;0 cos(a*t) sin(a*t);0 -sin(a*t) cos(a*t)];
%qaf=matlabFunction(qa);
%dqa=diff(qa,'t');
% daqf=matlabFunction(dqa);
% qad=[0 0 0;0 -a*sin(a*t) a*cos(a*t);0 -a*cos(a*t) -a*sin(a*t)];
qb=[cos(b*t) 0 -sin(b*t);0 1 0; sin(b*t) 0 cos(b*t)];
qc=[cos(c*t) sin(c*t) 0; -sin(c*t) cos(c*t) 0; 0 0 1];
q=qa*qb*qc;
dq=diff(q,'t');

dqf=matlabFunction(dq);

qf=matlabFunction(q);

qfT=(qf(tt))';
dqfT=(dqf(tt))';

bvc=[1.0; t^2; 0.5];
bvcf=matlabFunction(bvc);
dbv=diff(bvc,'t');
dbvf=matlabFunction(dbv);
bvcf1=bvcf(tt);
dbvf1=dbvf(tt);

% c1 = [5;3;4];
% t = 3;
% c2 = [101;30;2];
%%

xyz = reshape(xyz,65,65,129,3);
xyz = permute(xyz,[4 1 2 3]);

%%
xyz_new = zeros(3,65,65,129);

%%
%xyz_new = zeros(3,129,129,257);
for k=1:129
    for j=1:65
        for i=1:65
            
xyz_new(:,i,j,k) = qfT*xyz(:,i,j,k)-qfT*bvcf1;
        end
    end
end
%%
xyz_new = permute(xyz_new,[2 3 4 1]);
xyz_new = reshape(xyz_new,65*65*129*3,1);
fid = fopen('sullivan_newobj.xyz','w');
fwrite(fid,12,'integer*4');
fwrite(fid,65,'integer*4');
fwrite(fid,65,'integer*4');
fwrite(fid,129,'integer*4');
fwrite(fid,12,'integer*4');
fwrite(fid,65*65*129*3*8,'integer*4');
fwrite(fid,xyz_new,'real*8');
fwrite(fid,65*65*129*3*8,'integer*4');
fclose(fid);
%%
fid = fopen('sullivan.q','r');
fread(fid,6,'integer*4');
num = fread(fid,1,'integer*4');
sol = fread(fid,65*65*129*3,'real*8');
fclose(fid);
sol = reshape(sol,65,65,129,3);
sol = permute(sol,[4 1 2 3]);
sol_new = zeros(3,65,65,129);
for k=1:129
    for j=1:65
        for i=1:65
            c1=dqfT*xyz(:,i,j,k)-dqfT*bvcf1-qfT*dbvf1;
            
sol_new(:,i,j,k) = qfT*sol(:,i,j,k)+c1;
        end
    end
end

sol_new = permute(sol_new,[2 3 4 1]);
sol_new = reshape(sol_new,65*65*129*3,1);
fid = fopen('sullivan_newobj.q','w');
fwrite(fid,16,'integer*4');
fwrite(fid,65,'integer*4');
fwrite(fid,65,'integer*4');
fwrite(fid,129,'integer*4');
fwrite(fid,3,'integer*4');
fwrite(fid,16,'integer*4');
fwrite(fid,65*65*129*3*8,'integer*4');
fwrite(fid,sol_new,'real*8');
fwrite(fid,65*65*129*3*8,'integer*4');
fclose(fid);
%% 
%seeding points
p1=[0.02;0.05;2.0];p2=[0;-0.45;0.5];
p3=[0.45;0;0.5];p4=[-0.45;0;0.5];
temp=0.45*sqrt(2)/2;
p5=[temp;temp;0.5];p6=[-temp;-temp;0.5];
p7=[-temp;temp;0.5];p8=[temp;-temp;0.5];
%qfT*xyz(:,i,j,k)-qfT*bvcf1;
p88_new=[2.25; 1.41; 2.50];
p88=qfT'*p88_new+bvcf1;
p1_new = qfT*p1-qfT*bvcf1;
p2_new = qfT*p2-qfT*bvcf1;
p3_new =  qfT*p3-qfT*bvcf1;
p4_new = qfT*p4-qfT*bvcf1;
p5_new = qfT*p5-qfT*bvcf1;
p6_new = qfT*p6-qfT*bvcf1;
p7_new =  qfT*p7-qfT*bvcf1;
p8_new =  qfT*p8-qfT*bvcf1;

