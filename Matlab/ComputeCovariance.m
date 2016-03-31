%% Configuration
addpath(horzcat(pwd,'\Helpers\'));
param.cov=1;
param.diagrem=1;
param.embedd=1;
chromosomes=23;

%%
sizeG=zeros(1,chromosomes);
G=cell(chromosomes);
for chrid1 = 1:chromosomes
    for chrid2 = chrid1:chromosomes
        chr1=num2str(chrid1);
        chr2=num2str(chrid2);
        contigsLeft=load(horzcat('contigMap.',chr1,'.csv'));
        contigsRight=load(horzcat('contigMap.',chr2,'.csv'));
        n1=length(contigsLeft); n2=length(contigsRight);
        sizeG(chrid1)=n1; sizeG(chrid2)=n2;
        variantMatrix=load(horzcat('phaseMat_',param.resolution,'.',chr1,'.',chr2,'.csv'));
        newvals=variantMatrix(:,3);
        M=sparse(variantMatrix(:,1),variantMatrix(:,2),newvals,2*length(contigsLeft),2*length(contigsRight));
        if (isequal(chrid1,chrid2))
            cM=mat2cell(M, [n1 n1], [n2 n2]);
            t=double(param.diagrem);
            %diagonal preconditioning #1
            M=[triu(cM{1,1},t)+triu(cM{1,1},t)'+diag(diag(cM{1,1})) triu(cM{1,2},t)+triu(cM{2,1},t)';...
                (triu(cM{1,2},t)+triu(cM{2,1},t)')' triu(cM{2,2},t)+triu(cM{2,2},t)'+diag(diag(cM{2,2}))];
        end
        G{chrid1,chrid2}=M;
    end
end
% complete tril and low pass filter
for chrid1 = 2:chromosomes
    for chrid2 = 1:chrid1
        G{chrid1,chrid2}=G{chrid2,chrid1}';
    end
end
tG=cell2mat(G);
valcap=max(prctile(diag(tG),80),prctile(tG(tG>0),80));
tG(tG>valcap)=valcap;
tG=spdiags(repmat(valcap,length(tG),1),0,tG);
tG=mat2cell(tG,2*sizeG,2*sizeG);

% serialize covariance embeddings
for chrid1 = 1:chromosomes
    if(param.cov==true)
        x=[tG{chrid1,:}]';
        y=vertcat(tG{:,chrid1});
        M=x'*y;
        M=M-min(M(:));
    else
        M=tG{chrid1,chrid1};
    end
    filenm=horzcat('ChrCovEmbedd.',num2str(chrid1),'.mat');
    save(horzcat('tmp.',filenm),'M');
end
