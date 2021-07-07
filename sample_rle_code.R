#commonCols = sample(intersect(intersect(colnames(ruv)[-1],colnames(raw)),colnames(nsolver)),300)
commonCols = intersect(intersect(colnames(ruv)[-1],colnames(raw)),colnames(nsolver))
##### RLE PLOTS
rawCommon = raw[,commonCols,with=F]
rawCommon = log2(rawCommon[,-1]+1)
ruvCommon = ruv[,commonCols,with=F]
ruvCommon = ruvCommon[,-1]
nsolverCommon = nsolver[,commonCols,with=F]
nsolverCommon = nsolverCommon[,-1]

rawCommon = rawCommon - apply(rawCommon,1,median)
ruvCommon = ruvCommon - apply(ruvCommon,1,median)
nsolverCommon = nsolverCommon - apply(nsolverCommon,1,median)

stackRaw = stack(rawCommon)
stackRaw$Method = 'Raw'
stackRUV = stack(ruvCommon)
stackRUV$Method = 'RUVSeq'
stackN = stack(nsolverCommon)
stackN$Method = 'nSolver'

stack = rbind(stackRaw,stackRUV,stackN)
er = ruv.pData[,c('Sample','Phase123')]
colnames(stack)[2] = 'Sample'
stack = merge(stack,er,by='Sample')
stack$Phase = as.factor(stack$Phase123)
levels(stack$Phase) = c('P1 (1993-1996)',
                        'P2 (1996-2001)',
                        'P3 (2008-2013)')
cols = c(sample(unique(as.character(stack$Sample[stack$Phase == 'P1 (1993-1996)'])),100),
         sample(unique(as.character(stack$Sample[stack$Phase == 'P2 (1996-2001)'])),100),
         sample(unique(as.character(stack$Sample[stack$Phase == 'P3 (2008-2013)'])),100))
stack$Method = factor(stack$Method,
                      levels = c('Raw','nSolver','RUVSeq'))
stack = stack[order(stack$Phase),]
stack = subset(stack,Sample %in% cols)
stack$Sample = factor(stack$Sample,
                      levels = unique(stack$Sample))
rle_plots = ggplot(data = stack,aes(x = Sample,y = values, color = Phase)) +
    geom_boxplot(coef = 6) + facet_wrap(~Method,ncol = 1) + theme_minimal() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=24),
          plot.title = element_text(size = 30),
          legend.title=element_text(size=20),
          legend.text=element_text(size=20),
          strip.text = element_text(size=24),
          panel.spacing=unit(1, "lines"),
          panel.border = element_rect(color = "grey",
                                      fill = NA, size = .1),
          legend.position = 'bottom',
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + xlab('Sample') +
    ylab('Median deviation of log expression') + 
    scale_color_manual(values = c('forestgreen','gold4','purple')) +
    ylim(c(-4,4)) + labs(color = 'CBCS phase')