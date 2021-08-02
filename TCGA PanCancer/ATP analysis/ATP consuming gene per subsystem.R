###barplot ### x axis is subsystem, y is the number of gene (2021.03.17)
load(gene_to_HMRdata_df)
gene_to_subsystem_consuming = select(gene_to_HMRdata_df[gene_to_HMRdata_df$factor == 1, ], gene, SUBSYSTEM)
gene_to_subsystem_consuming = group_by(gene_to_subsystem_consuming, SUBSYSTEM)
gene_to_subsystem_consuming = summarise(gene_to_subsystem_consuming, num_gene = n())


# ggplot(gene_to_subsystem_consuming, aes(x=SUBSYSTEM, y=num_gene)) +
#     geom_bar(stat="identity", fill="steelblue", width=0.3) +
#     geom_text(aes(label=num_gene), vjust=-0.6, color="black", size=3) +
# #    geom_jitter(width=0.15) +
# #    ylim(0,40) +
#     #coord_cartesian((ylim))
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle=90,hjust=1, vjust=0.5)) +
#     xlab('the number of gene') +
#     ylab('Subsystem')

ggplot(gene_to_subsystem_consuming, aes(x=num_gene, y=SUBSYSTEM)) +
    geom_bar(stat="identity", fill="steelblue", width=0.3) +
    geom_text(aes(label=num_gene), vjust=0, hjust=-0.21, color="black", size=3) +
    #    geom_jitter(width=0.15) +
    #    ylim(0,40) +
    #coord_cartesian((ylim))
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90,hjust=0.5, vjust=0.5)) +
    xlab('Subsystem') +
    ylab('the number of ATP consuming gene')

save(gene_to_HMRdata_df, gene_to_RXNID_result, gene_to_subsystem_consuming, file='ATP_consuming_gene_per_subsystem.rda')


#filter : #gene >= 10
filter_subsystem_consuming =  gene_to_subsystem_consuming[gene_to_subsystem_consuming$num_gene>=10,]
