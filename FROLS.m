clear all
N = 200;
u = -1+2*rand(N,1);
y = zeros(N,1);
e = 0.1*randn(N,1);
%ny = 2;
%nu = 2;
%ne = 0;
l = 3;
del = 0.05;
y(1:2) = e(1:2);

for k = 3:N
    y(k) = -0.605*y(k-1)-0.163*y(k-2)^2+0.588*u(k-1)-0.24*u(k-2)+e(k);
end

z = 1;
for f = 0:l
    for h = 0:l
        for i = 0:l
            for j = 0:l
                if f+h+i+j>l
                    continue;
                end
                for k = 3:N
                    t(k-2) = (y(k-1)^j)*(y(k-2)^i)*(u(k-1)^h)*(u(k-2)^f);
                end
                %fprintf('z=%d y1=%d y2=%d u1=%d u2=%d \n',z,j,i,h,f);
                Z(z,:) = [j i h f];
                D(:,z) = t; %Dictionary of candidate model terms
                z = z+1;
            end
        end
    end
end

M = z-1;
y1 = y(3:N);
sig = y1'*y1;

%Step-1: s=1
for j = 1:M
    q = D(:,j);
    g1(j) = (y1'*q)/(q'*q);
    ERR(j) = g1(j)^2*(q'*q)/sig;
end
[c, b(1)] = max(ERR);
Al(:,1) = D(:,b(1)); %Matrix consisting of alpha
Q(:,1) = D(:,b(1)); %Matrix consisting of q
g(1) = g1(b(1)); %Vector consisting of model parameters, g
A(1,1) = 1;
err(1) = ERR(b(1)); %Vector consisting of ERR values
ers = err(1);

%Step-2: s>=2
for i = 2:M
    k = 1;
    x = zeros(M-i+1,1);
    g1 = zeros(1,M-i+1);
    D0 = zeros(N-2,M-i+1);
    Q1 = zeros(N-2,M-i+1);
    ERR = zeros(M-i+1,1);
    
    for j = 1:M
        if find(b==j)
            continue;
        else
            x(k) = j;
            k = k+1;
        end
    end
    
    k = 1;
    for j = 1:size(x)
        D0(:,k) = D(:,x(j));
        s = zeros(N-2,1);
        p = D(:,x(j));
        for r = 1:i-1
            q1 = Q(:,r);
            s = s+((p'*q1)/(q1'*q1)*q1);
        end
        q = p-s;
        g1(j) = (y1'*q)/(q'*q);
        ERR(j) = g1(j)^2*(q'*q)/sig;
        Q1(:,j) = q;
        k = k+1;
    end
    [c, d] = max(ERR);
    b(i) = x(d);
    Al(:,i) = D0(:,d); %Matrix consisting of alpha
    Q(:,i) = Q1(:,d); %Matrix consisting of q
    g(i) = g1(d); %Vector consisting of model parameters, g
    for r = 1:i-1
        A(r,i) = (Q(:,r)'*Al(:,i))/(Q(:,r)'*Q(:,r));
    end
    A(i,i) = 1;
    err(i) = ERR(d); %Vector consisting of ERR values
    ers = ers+err(i); %Sum of ERR values
    if (1-ers)<=del
        break;
    end
end

fprintf('Model terms:\n');
for k = 1:size(g,2)
    fprintf('y(k-1)^%d*y(k-2)^%d*u(k-1)^%d*u(k-2)^%d\n',Z(b(k),:));
end
fprintf('\n\nModel parameters:\n');
for k = 1:size(g,2)
    fprintf('%.4f\n',g(k));
end
