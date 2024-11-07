src_new = residual;
n = length(src_new);
Bi = [1/2 1/2];
Ai = 1;
Bo = [3/2 1];
m = 4;
tr_mat = zeros(n,m);

load('Abart.txt');
load('bird.txt')
load('ae-bat.txt')
S = Abart';
% S = bird';
% S = ae_bat';
N = length(S);
k = (S(1:end-1) - S(2:end)) ./ (S(1:end-1) + S(2:end));
ins = zeros(2,N);
outs = zeros(2,N);
uend = 0;
lend = 0;
R0 = 0.4;
RL = -0.1;
y = zeros(1,n);
zf = 0;
zfy = 0;
for i = 1:n
    uend = outs(1,N);
    lend = outs(2,1);
    
    outs(1,2:N) = ins(1,1:N-1).*(1+k) - k.*ins(2,2:N);
    outs(2,1:N-1) = ins(2,2:N).*(1-k) + k.*ins(1,1:N-1);
    outs(1,1) = lend*R0 + src_new(i);
%     outs(2,N) = uend*RL;
    [z,zf] = filter(Bi,Ai,uend,zf);
    outs(2,N) = -z;
    ins = outs;
    [y(i),zfy] = filter(Bo,Ai,uend,zfy);
    y(i) = uend;
end
figure
plot(y)
% tr_mat(:,1) = y';
% 
% % S = Abart';
% S = bird';
% % S = ae_bat';
% N = length(S);
% k = (S(1:end-1) - S(2:end)) ./ (S(1:end-1) + S(2:end));
% ins = zeros(2,N);
% outs = zeros(2,N);
% uend = 0;
% lend = 0;
% R0 = 0.4;
% RL = -0.1;
% y = zeros(1,n);
% zf = 0;
% zfy = 0;
% for i = 1:n
%     uend = outs(1,N);
%     lend = outs(2,1);
%     
%     outs(1,2:N) = ins(1,1:N-1).*(1+k) - k.*ins(2,2:N);
%     outs(2,1:N-1) = ins(2,2:N).*(1-k) + k.*ins(1,1:N-1);
%     outs(1,1) = lend*R0 + src_new(i);
% %     outs(2,N) = uend*RL;
%     [z,zf] = filter(Bi,Ai,uend,zf);
%     outs(2,N) = -z;
%     ins = outs;
%     [y(i),zfy] = filter(Bo,Ai,uend,zfy);
%     y(i) = uend;
% end
% tr_mat(:,end) = y';
% 
% for i = 1:n
%     l = linspace(tr_mat(i,1),tr_mat(i,end),m);
%     tr_mat(i,:) = l;
% end
% 
% y = reshape(tr_mat,1,[]);