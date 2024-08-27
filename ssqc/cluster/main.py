from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.cluster import KMeans, HDBSCAN
from sklearn.mixture import GaussianMixture


def cluster():
    ''' merge and cluster all QC report under input_dir into one file '''

    import argparse
    parser = argparse.ArgumentParser(description='Merge and cluster all QC reports in given directory')
    parser.add_argument('input_dir', type=str, help='Input directory')
    parser.add_argument('--output', type=str, default = "QC_Report.tsv", help='Output Filename')
    parser.add_argument('--k', type=int, default = 3, help='Number of Kmeans clusters')
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output = Path(args.output)

    result = pd.DataFrame()

    for file in input_dir.iterdir():
        if file.is_file() and file.suffix == ".tsv":
            df = pd.read_csv(file, sep="\t")
            result = pd.concat([result, df])

    
    # Extracting feature names by excluding the first column (which is filename, not a feature)
    #feature_names = df.columns[1:]
    feature_names = ['background','spikiness']
    print(f"Clustering based on {feature_names}")

    # Normalize the features
    scaler = RobustScaler()  # or use MinMaxScaler() if you prefer
    normalized_features = scaler.fit_transform(result[feature_names])
    print(normalized_features)

    # Cluster with Kmeans
    kmeans = KMeans(n_clusters=args.k, n_init='auto')
    result['kmeans'] = kmeans.fit_predict(normalized_features)

    # Cluster with Gaussian Mixture
    gmm = GaussianMixture(n_components=args.k)
    result['gaussian mixture'] = gmm.fit_predict(normalized_features)

    # Cluster with HDBSCAN
    dbscan = HDBSCAN(min_samples=5)
    result['HDBSCAN'] = dbscan.fit_predict(normalized_features)

    result.to_csv(output, sep="\t", index=False)
