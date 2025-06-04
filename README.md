# Spectral Compression Algorithms

This project implements several spectral algorithms for computing dominant eigenpairs of symmetric matrices, and applies them to image compression tasks in grayscale and color.

## 🔍 Project Structure

- `methods/`: Contains all algorithm implementations (e.g., power method, subspace iteration, etc.)
- `tests/`: Unit tests and comparison scripts
- `image_reconstruction/`: Scripts and data for image compression
- `reports/`: Final reports in French and English
- `figures/`: RMSE plots and reconstructed images

## 📊 Implemented Methods

- Power method (`v11`, `v12`)
- Subspace Iteration (`v0`, `v1`, `v2`, `v3`)
- Rayleigh-Ritz projection
- Comparison with MATLAB’s native `eig`

## 📸 Application

The methods are applied to compress grayscale and RGB images using truncated eigenvalue decomposition.

## 🛠️ Requirements

- MATLAB (tested with R2023a)

## 📄 Reports

The French version of the report is available in `reports/Rapport.pdf`, and the English version can be found in `Report.pdf`.

## 📈 Author

Ralph Khairallah  
ENSEEIHT 2024–2025

