function ensamble=ensambleEmbedding(file,isPaired,initialization)
    if(isstr(file) && exist(file,'file'))
        display(horzcat('loading data from',file));
        load(file);
    else
        M=file;
    end
    if(~exist('isPaired','var'))
        isPaired=false;
    end
    
    opts=statset('Display','iter','MaxIter',100,'TolX',10^-10,'tolfun',10^-10);%300);
    n=size(M,1);
    if (n>10000)
        cM=mat2cell(M,[n/2 n/2], [n/2 n/2]);
        first=[1:floor(n/4) ((n/2)+1):((n/2)+floor(n/4))];
        if(mod(n/2,2)==1)
            second=[ceil(n/4):(n/2) (n/2+ceil(n/4)):n];
        else
            second=[1+ceil(n/4):(n/2) (n/2+ceil(n/4)+1):n];
        end
        display('Split ensamble. Computing first.')
        
        if(exist('initialization','var'))
            ens1=ensambleEmbedding(M(first,first),isPaired,initialization(first,:));
            display('Computing second.')
            ens2=ensambleEmbedding(M(second,second),isPaired,initialization(second,:));
        else
            tM=M./max(M(:));
            tM=spdiags(ones(n,1),0,tM);
            initAll=cmdscale(tM,3);
            ens1=ensambleEmbedding(M(first,first),isPaired,initAll(first,:));
            display('Computing second.')
            ens2=ensambleEmbedding(M(second,second),isPaired,initAll(second,:));
        end
        display('Combining.')
        ensamble.coords=[cellfun(@(x) first(x),ens1.coords,'uniformoutput',false) cellfun(@(x) second(x),ens2.coords,'uniformoutput',false)];
        ensamble.M=[ens1.M ens2.M];
        ensamble.eigenvalues=[ens1.eigenvalues ens2.eigenvalues];
        ensamble.cmds=[ens1.cmds ens2.cmds];
        ensamble.mds=[ens1.mds ens2.mds];
        ensamble.stress=[ens1.stress ens2.stress];
        ensamble.embeded=[ens1.embeded ens2.embeded];
        ensamble.eigenvalues=[ens1.eigenvalues ens2.eigenvalues];
        ensamble.embeddingType=[ens1.embeddingType ens2.embeddingType];

    else
        M=M./max(M(:));
        M=spdiags(ones(n,1),0,M);

        blocks=components(M);
        [blockcount]=hist(blocks,unique(blocks));
        ensamble.coords=cell(1,max(blocks));
        ensamble.M=cell(1,max(blocks)); 
        ensamble.eigenvalues=cell(1,max(blocks));
        ensamble.cmds=cell(1,max(blocks));
        ensamble.mds=cell(1,max(blocks));
        ensamble.stress=zeros(1,max(blocks));
        ensamble.embeded=false(1,max(blocks));
        display(horzcat('found ',num2str(max(blocks))));
        for comp=1:max(blocks)
            boolids=blocks==comp;
            if(isPaired)
                %orids=xor(boolids(1:n/2) , boolids(((n/2)+1):n));
                %should we set a disconnected paired coord to epsilon? or keep as unknown?
                tbids=(boolids(1:n/2) | boolids(((n/2)+1):n));
                boolids=[tbids tbids];
            end
            ensamble.coords{comp}= find(boolids);
            ensamble.M{comp}=M(boolids,boolids);
            if(sum(boolids)>3) %embedding is relevant
                display(horzcat('embedding component ',num2str(comp), ' of size ',num2str(sum(boolids))));
                ensamble.embeded(comp)=true;
                if (exist('initialization','var'))
                    ensamble.cmds{comp}=initialization(boolids,:);
                else
                    [ensamble.cmds{comp},ensamble.eigenvalues{comp}]=cmdscale(M(boolids,boolids),3);
                end
                mat=full(M(boolids,boolids));
                mat(mat==0)=nan;
                try
                    try
                        [ensamble.mds{comp},ensamble.stress(comp)]=mdscale(ensamble.M{comp},3,...
                            'start',ensamble.cmds{comp}(:,1:3),'options',opts);
                        ensamble.embeddingType{comp}='cmdbased';
                    catch
                        [Y, stress] =mdscale(ensamble.M{comp},3,'options',opts,'start','random');
                        ensamble.embeddingType{comp}='random';
                    end
                catch
                    ensamble.embeddingType{comp}='failed';
                end
            else
                ensamble.embeddingType{comp}='none';
            end
        end
    end
    if(isstr(file) && exist(file,'file'))
        splt=strsplit(file,'.');
        save(strjoin(splt(2:end),'.'),'ensamble');
    end
end