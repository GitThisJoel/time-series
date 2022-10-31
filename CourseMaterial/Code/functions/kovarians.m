function [r,tau]=kovarians(C,A,m,method)
%[r,tau]=kovarians(C,A,m)
%
%  C=[1 c1 ... cq]
%  A=[1 a1 ... ap]
%  m is the maximum tau value.

% Finn Lindgren 2000-02-15 Preallokering av matriser,
%                          hastighetsvinst ca faktor 5.
%               2000-02-16 'filter' i.st.f. for-loop.

if nargin<4
  method='new';
end

C=C/A(1);
A=A/A(1);

p=length(A)-1;
q=length(C)-1;
n=max([p q]);
switch method
 case 'old'
  A1temp=[A zeros(1,n-p)]';
  C1temp=[C zeros(1,n-q)]';
  A2temp=A1temp;
  A1=A1temp;
  A2=zeros(n+1,1);
  for i=2:n+1
    A1temp=[0;A1temp];
    A2temp=[A2temp;0];
    A1=[A1 A1temp(1:n+1,1)];
    A2=[A2 A2temp(i:n+i,1)];
  end
  A3=A1(1:q+1,1:q+1);
  C1=C1temp;
  for i=2:q+1
    C1temp=[C1temp;0];
    C1=[C1 C1temp(i:n+i,1)];
  end
 case 'new'
  Atemp=[zeros(1,n) A zeros(1,2*n-p)]';
  A1=zeros(n+1,n+1);
  for i=1:n+1
    A1(:,i)=Atemp(n+2-i:n+2-i+n,1);
  end
  A2=zeros(n+1,n+1);
  for i=2:p+1
    A2(:,i)=Atemp(n+i:n+i+n,1);
  end
  A3=A1(1:q+1,1:q+1);
  Ctemp=[zeros(1,n) C zeros(1,2*n-q)]';
  C1=zeros(n+1,q+1);
  for i=1:q+1
    C1(:,i)=Ctemp(n+i:n+i+n,1);
  end
end

r=((A1+A2)\(C1*(A3\C')))';

if p>0
  switch method
   case 'old'
    for i=n+2:m+1
      r=[r -r(1,i-1:-1:i-p)*A(1,2:p+1)'];
    end
    r=r(1,1:m+1);
   case 'new'
    Atemp = [A zeros(1,n-p+1)];
    zi=filter(Atemp,1,r);
    r=filter(1,Atemp,zeros(1,m+1),zi);
  end
else
  r=[r zeros(1,m-length(r)+1)];
  r=r(1,1:m+1);
end

tau=0:m;
