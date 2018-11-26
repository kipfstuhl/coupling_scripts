r = 0.5;
h1 = 0.1;
l = 4;
p5 = [l/2, r + h1];
h2 = 0.5;
p6 = [l/2+l/8, r + 3*h1/2];
p7 = [13/16*l, r + 3*h1/2+h2/2];
p8 = [l, r + 3*h1/2+h2];
nnorm = sqrt( (3/16*l)*(3/16*l) + (h2/2)*(h2/2));
n1 = (3/16*l)/nnorm;
n2 = (h2/2)/nnorm;
r1 = 0.5;
t1 = n2;
t2 = -n1;
p9 = [l + t1*r1, r + 3*h1/2+h2 + t2*r1];
p10 = [13/16*l + t1*r1, r + 3*h1/2+h2/2 + t2*r1];
p11 = [l/2+l/8 + t1*r1, r + 3*h1/2 + t2*r1];
h3 = 0.00;
h4 = 0.1;
nnorm = sqrt( (3/16*l)*(3/16*l) + (h4/2)*(h4/2));
n1 = (3/16*l)/nnorm;
n2 = -(h4/2)/nnorm;
r2 = 0.3;
t1_ = -n2;
t2_ = +n1;
p12 = [(l/2+l/8 + t1*r1 + l/2+l/8+ t1_*r2)/2, (r + 3*h1/2 + t2*r1 -(r + 3*h3/2)+ t2_*r2)/2];

ni = 4;
nj = 3;
coefs = zeros(2,ni,nj);

coefs(:,1,1) = p12;
coefs(:,2,1) = p11;
coefs(:,3,1) = p10;
coefs(:,4,1) = p9;

coefs(:,1,2) = [1.75 0];
coefs(:,2,2) = (p11 + p6)/2;
coefs(:,3,2) = (p10 + p7)/2;
coefs(:,4,2) = (p9 + p8)/2;

coefs(:,1,3) = p5;
coefs(:,2,3) = p6;
coefs(:,3,3) = p7;
coefs(:,4,3) = p8;

nurb = nrbmak(coefs,{[0 0 0 0 1 1 1 1],[0 0 0.5 1 1]});
% nrbplot(nurb,[30 30],'colormap','white')
% view(0,90)
% hold on
% 
% for i = 1:ni
%     for j = 1:nj
%         p = coefs(:,i,j);
%         plot(p(1),p(2),'.r','Markersize',10);
%     end
% end