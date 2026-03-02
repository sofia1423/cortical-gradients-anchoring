% ================================================================================
% @file    dimensionality_reduction.m
% @brief   Dimensionality reduction and classification of structure data
% @author  [Yuqing Li]
% @date    2026-03-03
% @version 1.0
% @license MIT (See LICENSE file for details)
% @contact [202321250018@mail.bnu.edu.cn]
% @note    The 91-region structural data were embedded into 2d
% low-dimensional spaces using t-SNE, CMDS, or UMAP, then classified into
% different groups by k-means or GMM.


%% dimensionality reduction

clear all

% 1. Transfer sequence

cor_base=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];

% 2. generate dimensionality reduction maps

    % 2.1  t-SNE
    clear All_tsne;
    combname(1).name='tsne-kmeans';
    combname(2).name='tsne-GMM';
    combname(3).name='CMDS-kmeans';
    combname(4).name='CMDS-GMM';
    combname(5).name='UMAP-kmeans';
    combname(6).name='UMAP-GMM';


    kwhich=cor_base;
    load('cor.mat')
    cor_refine=cor_new;
    clear A;
    A=[1:1:size(cor_refine,1)];
    A_left = A(~ismember(A, kwhich));
    A_now=A(~ismember(A,A_left));
    clear cor;
    cor=cor_new;
    for i2=size(cor_new,1):-1:1
        if(sum(ismember(A_left,i2))==1);
            cor(:,i2)=[];
            cor(i2,:)=[];
        end;
    end;

    clear DDist;DDist=sqrt(2*(1-cor));%sqrt(2 * (1 - A))

    for i=1:size(DDist,1)
        DDist(i,i)=0;
    end;

    nsize=size(cor,1);

    % 1.dimensionality reduction
    dimensionality_number=2;

    Iteration_times=statset('MaxIter',1000);
    clear total_Dist;total_Dist(1:size(DDist,1),1:size(DDist,1))=0;
    clear dist_variance1;n=0;
    clear sum_relative_distance;nn=0;
    for p=1:500;
        clear Y;
        Y=tsne(DDist(:,1:size(DDist,2)),'NumDimensions',dimensionality_number,'Options',Iteration_times);
    %     Y=cmdscale(DDist(:,1:size(DDist,2)),dimensionality_number);
        clear Dist
        for i=1:size(DDist,1);
            for ii=1:size(DDist,1);
                clear dist;dist=0;
                for j=1:dimensionality_number;
                    dist=dist+power((Y(i,j)-Y(ii,j)),2);
                end;
                Dist(i,ii)=sqrt(dist);
                total_Dist(i,ii)=total_Dist(i,ii)+Dist(i,ii);
            end;
        end;
        n=n+1;
        nn=nn+1;
        dist_variance1(:,n)=total_Dist(:,1)./n;
        sum_relative_distance(1,n)=sum(sum(total_Dist(:,:)))./nn;
    end;
    total_Dist=total_Dist./500;
    sum_relative_distance=sum_relative_distance./((nsize-1)*(nsize-1)-(nsize-1));
    All_tsne{1,1}=total_Dist;

    save All_tsne.mat All_tsne


    % 2.2  CMDS
    clear All_CMDS;

    kwhich=cor_base;
    load('cor.mat')
    cor_refine=cor_new;
    clear A;
    A=[1:1:size(cor_refine,1)];
    A_left = A(~ismember(A, kwhich));
    A_now=A(~ismember(A,A_left));
    clear cor;
    cor=cor_new;
    for i2=size(cor_new,1):-1:1
        if(sum(ismember(A_left,i2))==1);
            cor(:,i2)=[];
            cor(i2,:)=[];
        end;
    end;

    clear DDist;DDist=sqrt(2*(1-cor));%sqrt(2 * (1 - A))

    for i=1:size(DDist,1)
        DDist(i,i)=0;
    end;

    nsize=size(cor,1);

    % 1.dimensionality reduction
    dimensionality_number=2;

    Iteration_times=statset('MaxIter',1000);
    clear total_Dist;total_Dist(1:size(DDist,1),1:size(DDist,1))=0;
    clear dist_variance1;n=0;
    clear sum_relative_distance;nn=0;
    for p=1:500;
        clear Y;
%             Y=tsne(DDist(:,1:size(DDist,2)),'NumDimensions',dimensionality_number,'Options',Iteration_times);
        Y=cmdscale(DDist(:,1:size(DDist,2)),dimensionality_number);
        clear Dist
        for i=1:size(DDist,1);
            for ii=1:size(DDist,1);
                clear dist;dist=0;
                for j=1:dimensionality_number;
                    dist=dist+power((Y(i,j)-Y(ii,j)),2);
                end;
                Dist(i,ii)=sqrt(dist);
                total_Dist(i,ii)=total_Dist(i,ii)+Dist(i,ii);
            end;
        end;
        n=n+1;
        nn=nn+1;
        dist_variance1(:,n)=total_Dist(:,1)./n;
        sum_relative_distance(1,n)=sum(sum(total_Dist(:,:)))./nn;
    end;
    total_Dist=total_Dist./500;
    sum_relative_distance=sum_relative_distance./((nsize-1)*(nsize-1)-(nsize-1));
    All_CMDS{1,1}=total_Dist;
    end;

    save All_CMDS.mat All_CMDS


    % 2.3  first part of UMAP
    clear DDist_UMAP;

    kwhich=cor_base;
    load('cor.mat')
    cor_refine=cor_new;
    clear A;
    A=[1:1:size(cor_refine,1)];
    A_left = A(~ismember(A, kwhich));
    A_now=A(~ismember(A,A_left));
    clear cor;
    cor=cor_new;
    for i2=size(cor_new,1):-1:1
        if(sum(ismember(A_left,i2))==1);
            cor(:,i2)=[];
            cor(i2,:)=[];
        end;
    end;

    clear DDist;DDist=sqrt(2*(1-cor));%sqrt(2 * (1 - A))

    for i=1:size(DDist,1)
        DDist(i,i)=0;
    end;

    DDist_UMAP{1,1}=DDist;

    save DDist_UMA.mat DDist_UMAP
    % Transfer to python code: code_for_UMAP.ipynb


%% 2. classification

clear all
cor_base=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];;


    % 2.1  k_means

    clear All_samegroup_matrix;

    kwhich=cor_bas;
    load('cor.mat')
    cor_refine=cor_new;
    clear A;
    A=[1:1:size(cor_refine,1)];
    A_left = A(~ismember(A, kwhich));
    A_now=A(~ismember(A,A_left));
    clear cor;
    cor=cor_new;
    for i2=size(cor_new,1):-1:1
        if(sum(ismember(A_left,i2))==1);
            cor(:,i2)=[];
            cor(i2,:)=[];
        end;
    end;

    clear DDist;DDist=sqrt(2*(1-cor));%sqrt(2 * (1 - A))

    for i=1:size(DDist,1)
        DDist(i,i)=0;
    end;

    nsize=size(cor,1);

    % 1.dimensionality reduction
    dimensionality_number=2;

    for i2=1:3;
        clear tcell;
        clear total_Dist;
        if i2==1;
            tcell=load('All_tsne.mat').All_tsne;
            total_Dist=tcell{1,1};
            for i3=1:2;
                if i3==1;
                    groupnumber=3;
                    kmeans_times=500;
                    clear Samegroup_matrix;Samegroup_matrix(1:size(DDist,1),1:size(DDist,1))=0;
                    for ii=1:kmeans_times;
                        clear idx;
                        [idx,C,sumd,D]=kmeans(total_Dist(:,1:size(DDist,1)),groupnumber);
                        for j=1:size(DDist,1);
                            k1=idx(j,1);
                            for jj=1:size(DDist,1);
                                k2=idx(jj,1);
                                if(k1==k2)
                                    Samegroup_matrix(j,jj)=Samegroup_matrix(j,jj)+1;
                                end;
                            end;
                        end;
                    end;
                    Samegroup_matrix=Samegroup_matrix./kmeans_times;
                    All_samegroup_matrix{1,(i2-1)*2+i3}=Samegroup_matrix;
                elseif i3==2;
                    ttotal_Dist=total_Dist;
                    ttotal_Dist(:,nsize)=[];
                    groupnumber=3;
                    kmeans_times=500;
                    clear Samegroup_matrix;Samegroup_matrix(1:size(DDist,1),1:size(DDist,1))=0;
                    for ii=1:kmeans_times;
                        clear idx;
                        gm = fitgmdist(ttotal_Dist, 3, 'RegularizationValue', 0.1);
                        idx = cluster(gm, ttotal_Dist);
                        for j=1:size(DDist,1);
                            k1=idx(j,1);
                            for jj=1:size(DDist,1);
                                k2=idx(jj,1);
                                if(k1==k2)
                                    Samegroup_matrix(j,jj)=Samegroup_matrix(j,jj)+1;
                                end;
                            end;
                        end;
                    end;
                    Samegroup_matrix=Samegroup_matrix./kmeans_times;
                    All_samegroup_matrix{1,(i2-1)*2+i3}=Samegroup_matrix;
                end;
            end;
        elseif i2==2;
            tcell=load('All_CMDS.mat').All_CMDS;
            total_Dist=tcell{1,1};
            for i3=1:2;
                if i3==1;
                    groupnumber=3;
                    kmeans_times=500;
                    clear Samegroup_matrix;Samegroup_matrix(1:size(DDist,1),1:size(DDist,1))=0;
                    for ii=1:kmeans_times;
                        clear idx;
                        [idx,C,sumd,D]=kmeans(total_Dist(:,1:size(DDist,1)),groupnumber);
                        for j=1:size(DDist,1);
                            k1=idx(j,1);
                            for jj=1:size(DDist,1);
                                k2=idx(jj,1);
                                if(k1==k2)
                                    Samegroup_matrix(j,jj)=Samegroup_matrix(j,jj)+1;
                                end;
                            end;
                        end;
                    end;
                    Samegroup_matrix=Samegroup_matrix./kmeans_times;
                    All_samegroup_matrix{1,(i2-1)*2+i3}=Samegroup_matrix;
                elseif i3==2;
                    ttotal_Dist=total_Dist;
                    ttotal_Dist(:,nsize)=[];
                    groupnumber=3;
                    kmeans_times=500;
                    clear Samegroup_matrix;Samegroup_matrix(1:size(DDist,1),1:size(DDist,1))=0;
                    for ii=1:kmeans_times;
                        clear idx;
                        gm = fitgmdist(ttotal_Dist, 3, 'RegularizationValue', 0.1);
                        idx = cluster(gm, ttotal_Dist);
                        for j=1:size(DDist,1);
                            k1=idx(j,1);
                            for jj=1:size(DDist,1);
                                k2=idx(jj,1);
                                if(k1==k2)
                                    Samegroup_matrix(j,jj)=Samegroup_matrix(j,jj)+1;
                                end;
                            end;
                        end;
                    end;
                    Samegroup_matrix=Samegroup_matrix./kmeans_times;
                    All_samegroup_matrix{1,(i2-1)*2+i3}=Samegroup_matrix;
                end;
            end;
        else 
            tcell=load('All_UMAP.mat').total_Dist_all';
            total_Dist=tcell{1,1};
            for i3=1:2;
                if i3==1;
                    groupnumber=3;
                    kmeans_times=500;
                    clear Samegroup_matrix;Samegroup_matrix(1:size(DDist,1),1:size(DDist,1))=0;
                    for ii=1:kmeans_times;
                        clear idx;
                        [idx,C,sumd,D]=kmeans(total_Dist(:,1:size(DDist,1)),groupnumber);
                        for j=1:size(DDist,1);
                            k1=idx(j,1);
                            for jj=1:size(DDist,1);
                                k2=idx(jj,1);
                                if(k1==k2)
                                    Samegroup_matrix(j,jj)=Samegroup_matrix(j,jj)+1;
                                end;
                            end;
                        end;
                    end;
                    Samegroup_matrix=Samegroup_matrix./kmeans_times;
                    All_samegroup_matrix{1,(i2-1)*2+i3}=Samegroup_matrix;
                elseif i3==2;
                    ttotal_Dist=total_Dist;
                    ttotal_Dist(:,nsize)=[];
                    groupnumber=3;
                    kmeans_times=500;
                    clear Samegroup_matrix;Samegroup_matrix(1:size(DDist,1),1:size(DDist,1))=0;
                    for ii=1:kmeans_times;
                        clear idx;
                        gm = fitgmdist(ttotal_Dist, 3, 'RegularizationValue', 0.1);
                        idx = cluster(gm, ttotal_Dist);
                        for j=1:size(DDist,1);
                            k1=idx(j,1);
                            for jj=1:size(DDist,1);
                                k2=idx(jj,1);
                                if(k1==k2)
                                    Samegroup_matrix(j,jj)=Samegroup_matrix(j,jj)+1;
                                end;
                            end;
                        end;
                    end;
                    Samegroup_matrix=Samegroup_matrix./kmeans_times;
                    All_samegroup_matrix{1,(i2-1)*2+i3}=Samegroup_matrix;
                end;
            end;
        end;
    end;

    save All_samegroup_matrix.mat All_samegroup_matrix


%% 3. Whether afferent is the new center

clear all

cor_base=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];

combname(1).name='tsne-kmeans';
combname(2).name='tsne-GMM';
combname(3).name='CMDS-kmeans';
combname(4).name='CMDS-GMM';
combname(5).name='UMAP-kmeans';
combname(6).name='UMAP-GMM';

load('All_samegroup_matrix.mat');

clear intragroup_corr;

kwhich=cor_bas;
load('cor.mat')
cor_refine=cor_new;
clear A;
A=[1:1:size(cor_refine,1)];
A_left = A(~ismember(A, kwhich));
A_now=A(~ismember(A,A_left));
clear cor;
cor=cor_new;
for i2=size(cor_new,1):-1:1
    if(sum(ismember(A_left,i2))==1);
        cor(:,i2)=[];
        cor(i2,:)=[];
    end;
end;

clear DDist;DDist=sqrt(2*(1-cor));%sqrt(2 * (1 - A))

for i=1:size(DDist,1)
    DDist(i,i)=0;
end;

nsize=size(cor,1);

dimensionality_number=2;

figure 
for i2=1:6;% 6 combination of dimensionality-reduction and classification

    clear Samegroup_matrix;
    Samegroup_matrix=All_samegroup_matrix{1,i2};
    clear Ci;
    Ci=modularity_dir(Samegroup_matrix);

    if(numel(unique(Ci))~=3)
        flist2(1,i2)=0;
    else

        for ttnum=1:size(A_now,2)
            kk2=A_now(1,ttnum);
            if(kk2==16)
                break;
            end;
        end;

        clear nCi;
        k1=Ci(4,1);
        k2=Ci(10,1);
        k3=Ci(ttnum,1);
        if(k1==k2)
            flist2(1,i2)=0;
        else
            for i=1:size(Ci,1)
                k=Ci(i,1);
                if(k==k1)
                    nCi(i,1)=1;
                elseif(k==k2)
                    nCi(i,1)=2;
                elseif(k==k3)
                    nCi(i,1)=3;
                else
                    nCi(i,1)=4;
                end;
            end;    

            clear Ci;
            Ci=nCi;
            nidx=Ci;

            nidx(:,2)=A_now;
            nidx(:,3)=[1:1:size(nidx,1)]';

            clear ridx;n=0;
            for i=1:max(nidx(:,1))
                for ii=1:size(nidx,1)
                    if(i==nidx(ii,1))
                        n=n+1;
                        ridx(n,:)=nidx(ii,:);
                    end;
                end;
            end;

            clear rridx;rridx(:,1:3)=0;n=0;
            for i=1:max(Ci);
                clear p;nn=0;
                for ii=1:size(ridx,1);
                    if(i==ridx(ii,1))
                        nn=nn+1;
                        p(nn,:)=ridx(ii,:);
                    end;
                end;

                clear pp;
                pp=p;

                if(i==1)
                    rridx(1:size(pp,1),:)=pp(1:size(pp,1),:);
                else
                    rridx((size(rridx,1)+1):(size(rridx,1)+size(pp,1)),:)=pp(1:size(pp,1),:);
                end;
            end;

            clear rcor_refine;
            for i=1:size(DDist,1);
                k1=rridx(i,2);
                for ii=1:size(DDist,1);
                    k2=rridx(ii,2);
                    rcor_refine(i,ii)=cor_refine(k1,k2);
                end;
            end;

            clear rSamegroup_matrix;
            for i=1:size(DDist,1);
                k1=rridx(i,3);
                for ii=1:size(DDist,1);
                    k2=rridx(ii,3);
                    rSamegroup_matrix(i,ii)=Samegroup_matrix(k1,k2);
                end;
            end;

            % correlation matrix
            subplot(3,2,i2);
            imagesc(rcor_refine)
            hold on

            for i=1:size(rridx,1)
                rridx(i,3)=i;
            end;

            clear rangeidx;
            for i=1:max(Ci);
                clear plist;p=0;
                for ii=1:size(rridx,1)
                    k=rridx(ii,1);
                    if(i==k);
                        p=p+1;
                        plist(p,:)=rridx(ii,:);
                    end;
                end;
                tmax=max(plist(:,3));
                tmin=min(plist(:,3));
                rangeidx(i,1)=tmin;
                rangeidx(i,2)=tmax;
            end;

            for i=1:max(Ci);
                k1=rangeidx(i,1);
                k2=rangeidx(i,2);
                rectangle('Position',[(k1-0.5) (k1-0.5) (k2-k1+1) (k2-k1+1)],'edgecolor','r','LineWidth',1);
            end;

            clear average_corr;average_corr(1:17,1:3)=0;
            clear I;
            for i=1:max(Ci);
                k=rangeidx(i,2)-rangeidx(i,1)+1;
                clear nlist;nlist(1:k,1)=0;
                for j=1:k;
                    for jj=1:k
                        nlist(j,1)=nlist(j,1)+rcor_refine(rangeidx(i,1)+j-1,rangeidx(i,1)+jj-1);
                    end;
                    nlist(j,2)=rridx(rangeidx(i,1)+j-1,2);
                end;
                for v=1:size(nlist,1)
                    average_corr(v,i)=nlist(v,1)/k;
                end;
                %the max value can be not only
                allmax=find(nlist==max(nlist,[],1));
                n=0;
                while(n<size(allmax,1));
                    b=allmax(n+1,1);
                    if(rridx(rangeidx(i,1)+b-1,1)==i)
                        I(i,1)=rridx(rangeidx(i,1)+b-1,2);
                        break;
                    end;
                    n=n+1;
                end;
            end;

            clear slist;
            for i=1:max(Ci)
                clear temp_corr;ntemp=0;
                for ii=1:size(rridx,1)
                    k1=rridx(ii,1);
                    s1=rridx(ii,2);
                    if(k1==i)
                        for iii=1:size(rridx,1)
                            if(ii~=iii)
                                k2=rridx(iii,1);
                                s2=rridx(iii,2);
                                if(k2==i)
                                    ntemp=ntemp+1;
                                    temp_corr(ntemp,1)=cor_refine(s1,s2);
                                    slist(ntemp,1)=s1;
                                    slist(ntemp,2)=s2;
                                end;
                            end;
                        end;
                    end;
                end;
                intragroup_corr{1,i2}(i,1)=sum(temp_corr)/ntemp;
            end;


            for i=1:size(I,1)
                k=I(i,1);
                for ii=1:size(rridx,1)
                    if(rridx(ii,2)==k)
                        kk=ii;
                    end;
                end;
                rectangle('Position',[(kk-0.5) (kk-0.5) 1 1],'edgecolor','r','LineWidth',1);
                text(kk+1,kk,num2str(k),'Color','red');
            end;


            colorbar
            xlim([1 size(DDist,1)])
            xticks([1:size(DDist,1)])
            % xticklabels(num2str(rridx(:,2)));

            xticklabels(num2str(rridx(:,2)));
            ylim([1 size(DDist,1)])
            yticks([1:size(DDist,1)])
            yticklabels(new_name(rridx(:,2)));

            hold off

            title(combname(i2).name, 'FontSize', 12, 'Color', 'red')  % 彩色标题

            % Here is for scatter plot in dimensionality reductd spcace
            
            new_idx= A_now';

            clear options;clear Location;clear k_rridx;
            options=statset('MaxIter',1000);
            Location=tsne(DDist(:,1:size(DDist,2)),'NumDimensions',dimensionality_number,'Options',options);
            k_rridx=kmeans(DDist(:,1:size(DDist,2)),3);
            % scatter figures
            % subplot(3,2,i2);
            % hold on
            % for i=1:size(Location,1)
            %     k=rridx(i,2);
            %     n=0;
            %     for ii=1:size(I,1);
            %         if(k==I(ii,1))
            %             n=n+1;
            %         end;
            %     end;
            %     if(n>0)
            %         scatter(Location(i,1),Location(i,2),120,k_rridx(i,1),'filled','d')
            %     else
            %         scatter(Location(i,1),Location(i,2),60,k_rridx(i,1),'filled','o')
            %     end;
            % end;
            % hold off
            % text(Location(:,1)+1,Location(:,2)+1,num2str(new_idx(:,1)));

        end;

    end;
end;
