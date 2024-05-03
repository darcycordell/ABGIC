function [B] = insert_arbitrary_RowsColumns(A,col)

M = size(A,1);
ins = length(col);
N = M+ins;
B = zeros(N);

indout = setxor(1:N,col);

B(indout,indout) = A;

% for i = 1:length(col)
%     A = [A(:,1:col(i)-1) zeros(size(A,1),1) A(:,col(i):end)];
%     A = [A(1:col(i)-1,:); zeros(1,size(A,2)); A(col(i):end,:)];
%     
% end

