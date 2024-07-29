
import pandas as pd
from pathlib import Path
import sys
import plotly.express as px
import pyranges as pr
import numpy as np
from sklearn.decomposition import PCA

class GenomeAnnot:
    genome_map = {'CP015399.2': 'YL32',
                  'CP015400.2': 'KB18',
                  'CP015401.2': 'I48',
                  'CP015402.2': 'YL27',
                  'CP015403.2': 'YL45',
                  'CP015404.2': 'I46',
                  'CP015405.2': 'YL58',
                  'CP015406.2': 'YL31',
                  'CP015407.2': 'YL2',
                  'CP015408.2': 'I49',
                  'CP015409.2': 'YL44',
                  'CP015410.2': 'KB1',
                  'GCF_000364265': 'ASF519',
                  'FQ312003.1': 'SL1344',
                  'FQ312003.1;FQ312003.1': 'SL1344',
                  'HE654725.1': 'SL1344',
                  'HE654726.1': 'SL1344',
                  'HE654724.1': 'SL1344',
                  'contig_15': 'contig_15',
                  'contig_21': 'contig_21',
                  'contig_26': 'contig_26',
                  'contig_46': 'contig_46',
                  'AQFU02000001.1': 'ASF 502',
                  'AQFU02000002.1': 'ASF 502',
                  'AQFU02000003.1': 'ASF 502',
                  'CP097573.1': 'ASF500',
                  'NZ_CP097810.1': 'ASF356',
                  'NZ_AQFR02000001.1': 'ASF360',
                  'NZ_AQFR02000002.1': 'ASF360',
                  'NZ_AQFR02000003.1': 'ASF360',
                  'NZ_CP097561.1': 'ASF361',
                  'NZ_AQFT02000001.1': 'ASF492',
                  'NZ_AQFT02000002.1': 'ASF492',
                  'NZ_AQFT02000003.1': 'ASF492',
                  'NZ_AQFT02000004.1': 'ASF492',
                  'NZ_AQFT02000005.1': 'ASF492',
                  'NZ_AQFT02000006.1': 'ASF492',
                  'NZ_AQFV02000001.1': 'ASF519',
                  'NZ_AQFV02000002.1': 'ASF519',
                  'NZ_AQFV02000003.1': 'ASF519',
                  'NZ_AQFV02000004.1': 'ASF519',
                  'NZ_AQFV02000005.1': 'ASF519',
                  'NZ_AQFV02000006.1': 'ASF519',
                  'NZ_CP097562.1': 'ASF457',
                  }
    annotation_columns = ['Chromosome',  'Feature', 'Start', 'End', 'Strand', 'ID',
                          'Name', 'locus_tag', 'gene_biotype', 'product']

    def __init__(self, gff_file, chr_names=pd.DataFrame()):
        self.gff_file = gff_file
        self.chr_names = chr_names
        self.feature = "gene"
        self.annot = self.process_gff()

    def process_gff(self):
        gff = pr.read_gff3(self.gff_file).as_df()[self.annotation_columns]
        return gff[gff['Feature'] == self.feature]

    def annotate_df(self, df):
        chr_map = pd.DataFrame.from_dict(self.genome_map, orient='index').reset_index()
        chr_map.columns = ['Chromosome', 'genome']
        if not self.chr_names.empty:
            chr_map = pd.concat([chr_map, self.chr_names])
        fdf = df.merge(self.annot, on='ID', how='left')
        fdf = fdf.merge(chr_map, on='Chromosome', how='left')
        return fdf


class CountDataSet:
    def __init__(self, data_dir):
        self.data_dir = Path(data_dir)
        self.count_data = pd.DataFrame()

    def load_count_files(self):
        pass


class DEResults:
    pass



# Currently not looking at these


class HtseqCounts(CountDataSet):
    count_col = 'count'

    def load_count_files(self):
        files = list(self.data_dir.rglob('*.txt'))
        df_list = []
        for f in files:
            df = pd.read_table(f, names=['Name', 'count'], header=None).assign(
                sample_id=f.stem.split(".")[0])
            df['Name'] = df.Name.str.split("gene-", expand=True)[1]
            df = df.dropna(subset=['Name'])
            df['genome'] = [self.genome_map.get(name.split(
                "_")[0], 'SL1344') for name in df.Name.values]
            df_list.append(df)
        self.count_data = (pd.concat(df_list).rename({self.count_col: 'read_counts',
                                                     'Name': 'locus_tag'}, axis=1)
                           .merge(self.annot, on='locus_tag', how='left'))


class SalmonCounts(CountDataSet):
    count_col = 'NumReads'
    gene_col = 'Name'

    def load_count_files(self):
        files = list(self.data_dir.rglob('quant.sf'))

        df_list = []
        for f in files:
            name = f.parent.stem.split("_quant")[0]
            print(name)
            df = pd.read_table(f).assign(sample_id=name)
            df = df.rename(
                {self.count_col: 'salmon_read_counts'}, axis=1)
            df['ID'] = (df[self.gene_col].str.split('ID=', expand=True)[1]
                        .str.split(";", expand=True)[0])
            df = df.drop(columns=[self.gene_col])
            df_list.append(df)
        self.count_data = pd.concat(df_list)


class FeatureCounts(CountDataSet):
    count_col = 'fc_read_counts'
    gene_col = 'ID'

    def load_count_files(self):
        files = list(self.data_dir.rglob("*.count.txt"))
        df_list = []
        print(files)
        for f in files:
            name = f.stem.split(".count")[0]
            print(name)
            df = pd.read_table(f, comment='#').assign(sample_id=name)
            df.columns = [self.gene_col, 'chr', 'start', 'end',
                          'strand', 'length', self.count_col, 'sample_id']
            df = df[['ID', 'length', self.count_col, 'sample_id']]
            df_list.append(df)
        self.count_data = pd.concat(df_list)

    def calculate_tpms(self, log=True):
        self.count_data['tpms'] = self.count_data[self.count_col]/self.count_data['length']*1000
        self.count_data['tpms'] = self.count_data['tpms']/self.count_data.groupby('sample_id')['tpms'].transform('sum')*1e6
        if log:
            self.count_data['tpms'] = np.log2(self.count_data['tpms'] + 0.5)
            
    @property
    def summary_df(self):
        """
        
        def read_featurcounts_summary(file):
            df = pd.read_table(file)
            df.columns = ['status', file.stem.split('.')[0]]
            df = df.set_index('status').T
            return df

        def get_featurecounts_summary(sum_dir):
            files = list(sum_dir.rglob("*.summary"))
            dfs = []
            for file in files:
                dfs.append(read_featurcounts_summary(file))
            return pd.concat(dfs)
        
        
        """
        files = list(self.data_dir.rglob("*.summary"))
        df_list = []
        for f in files:
            df = pd.read_table(f)
            name = df.columns[1].split("/")[-1].split('.')[0]
            df = df.assign(sample_id=name)
            df.columns = ['status', 'read_counts', 'sample_id']
            df_list.append(df)
        fdf = pd.concat(df_list)
        summary = fdf.groupby('sample_id').read_counts.sum().reset_index()
        summary.columns = ['sample_id', 'total']
        summary = (summary.merge(fdf[fdf.status == 'Assigned'][['read_counts', 'sample_id']], on='sample_id')
                   .rename({'read_counts': 'assigned'}, axis=1)
                   .merge(fdf[fdf.status == 'Unassigned_Unmapped'][['read_counts', 'sample_id']], on='sample_id')
                   .rename({'read_counts': 'unmapped'}, axis=1)
                   .merge(fdf[fdf.status == 'Unassigned_NoFeatures'][['read_counts', 'sample_id']], on='sample_id')
                   .rename({'read_counts': 'no_feature'}, axis=1)
                   .merge(fdf[fdf.status == 'Unassigned_Ambiguity'][['read_counts', 'sample_id']], on='sample_id')
                   .rename({'read_counts': 'ambiguous'}, axis=1))
        
        summary['percent_assigned'] = summary['assigned']/summary['total']*100
        summary['percent_unmapped'] = summary['unmapped']/summary['total']*100
        summary['percent_ambiguous'] = summary['ambiguous']/summary['total']*100
        summary['percent_no_feature'] = summary['no_feature'] / \
            summary['total']*100
        return summary


class SushiCounts(CountDataSet):
    count_col = "total_insertcount"
    gene_col = "#reference"

    def load_count_files(self):
        files = list(self.data_dir.rglob("*ushicounts"))
        df_list = []
        for f in files:
            name = f.stem.split(".")[0]
            print(name)
            df = pd.read_table(
                f, usecols=[0, 2, 6, 7, 8]).assign(sample_id=name)
            df = df.rename(columns={self.count_col: "sushi_insertcount"})
            df['ID'] = df['#reference'].str.split(
                ';', expand=True)[0].str.split('ID=', expand=True)[1]
            df = df.drop(columns=[self.gene_col])
            df_list.append(df)
        self.count_data = pd.concat(df_list)

        # self.count_data = fdf.merge(self.annot, on='ID', how='left')
        # self.count_data["genome"] = self.count_data['Chromosome'].replace(
        #     self.genome_map)



#PCA
        

def find_PCs(count_data, sampleData, num_pcs=2, num_genes=None, choose_by='variance'):
    """
    :param countData: each column is a sampleID, index is featureID
    :param sampleData:
    :param numPCs:
    :param numGenes:
    :return:
    """
    if num_genes:
        # calculate var for each, pick numGenes top var across samples -> df
        if choose_by == 'variance':
            genes = count_data.var(axis=1).sort_values(ascending=False).head(num_genes).index
            df = count_data.loc[genes].T
        else:
            pass
            # todo implement log2fc selection
    else:
        df = count_data.T
    pca = PCA(n_components=num_pcs)
    principalComponents = pca.fit_transform(df)
    pcs = [f'PC{i}' for i in range(1, num_pcs+1)]
    pDf = (pd.DataFrame(data=principalComponents, columns=pcs)
           .set_index(df.index))
    pc_var = {pcs[i]: round(pca.explained_variance_ratio_[i] * 100, 2) for i in range(0, num_pcs)}
    pDf2 = pDf.merge(sampleData, left_index=True, right_index=True)
    return pDf2, pc_var


#PLOTS
def composition_plot(df, x_col='sample_id', y_col='proportion', color_col='genome', 
                     genome_color_map={}, w=800, h=800):
    fig = px.bar(df, x=x_col, 
                 y=y_col, color=color_col, 
                 color_discrete_map=genome_color_map, 
                 template='plotly_white', 
                 width=w, height=h, 
                 hover_data=df.columns,
                 labels = {y_col: 'Transcriptome abundance', x_col:''})

    return fig