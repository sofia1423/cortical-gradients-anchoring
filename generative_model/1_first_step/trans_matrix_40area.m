function B=trans_matrix_40area(A,x,y)
% load('C:\Users\Ale\Nutstore\1\workfiles\work\macaque\SLN\SLN91.mat')
% SLN=SLN91';
% zero_columns = find(all(SLN == 0, 1));
load('Z:\work\spower_all\figure_7\generative_2026_119\MatF_SLN2021.mat')
zero_columns = find(all(MatF== 0, 1));

if x ==40;
    A(zero_columns,:)=[];
end
if y ==40;
    A(:,zero_columns)=[];
end
B=A;