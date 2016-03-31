for chrid1 = 1:chromosomes
    filenm=horzcat('ChrCovEmbedd.',num2str(chrid1),'.mat');
    system(horzcat('matlab -nosplash -nodesktop -r "ensambleEmbedding(''tmp.',...
        filenm,''',true); exit;" -logfile debug.',num2str(chrid1),'.log -sd ',pwd, ' &'));
end
