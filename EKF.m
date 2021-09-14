clear
N = 100;
x = [1;1];
xi = x;
P = 1000*eye(2);

q1 = 0.15*randn(N,1);
q2 = 0.15*randn(N,1);
r = 0.2*randn(N,1);
Q = [var(q1) 0; 0 var(q2)];
R = var(r);

%Actual states calculation
for i=1:N
    xa = [x(1)/(1+x(2)^2)+q1(i); x(1)*x(2)/(1+x(2)^2)+q2(i)];
    y(i) = xa(1)+r(i);
    xa1(i,:) = xa(1);
    xa2(i,:) = xa(2);
    x = xa;
end

%Estimated states calculation
for i=1:N
    A = [1/(1+xi(2)^2) -2*xi(1)*xi(2)/(1+xi(2)^2)^2;
         xi(2)/(1+xi(2)^2) xi(1)*(1-xi(2)^2)/(1+xi(2)^2)^2];
    x = [xi(1)/(1+xi(2)^2); xi(1)*xi(2)/(1+xi(2)^2)];
    P = A*P*A'+Q;
    C = [1 0]; 
    K = P*C'*pinv(C*P*C'+R);
    x = x+K*(y(i)-x(1));
    P = (eye(2)-K*C)*P;
    xk1(i,:) = x(1);
    xk2(i,:) = x(2);
    yk(i) = C*x;
    if sum((xa1(i)-xk1(i))^2+(xa2(i)-xk2(i))^2)/sum(xa1(i)^2+xa2(i)^2) < 0.01
        i1 = i;
        break;
    end
    xi = x;
end

t1 = 1:i1;
subplot(2,1,1)
plot(t1,xa1(1:i1),'m',t1,xk1,'b');
xlabel('k');
ylabel('x1');
legend({'Actual State','Estimated State'})
subplot(2,1,2)
plot(t1,xa2(1:i1),'m',t1,xk2,'b');
xlabel('k');
ylabel('x2');
legend({'Actual State','Estimated State'})
