#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="RAAVioliShort_env"
YML_FILE="raaviolishort_env.yml"

# ---------------- CREATE/UPDATE ENV ---------------- #
if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  echo "[INFO] Updating existing env '$ENV_NAME'..."
  conda env update -n "$ENV_NAME" -f "$YML_FILE" --prune
else
  echo "[INFO] Creating env '$ENV_NAME'..."
  conda env create -n "$ENV_NAME" -f "$YML_FILE"
fi

# ---------------- ACTIVATE ENV ---------------- #
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

# ---------------- DISABLE USER-SITE ---------------- #
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
mkdir -p "$CONDA_PREFIX/etc/conda/deactivate.d"

cat > "$CONDA_PREFIX/etc/conda/activate.d/disable_usersite.sh" <<'EOF'
export PYTHONNOUSERSITE=1
unset PYTHONPATH
EOF

cat > "$CONDA_PREFIX/etc/conda/deactivate.d/disable_usersite.sh" <<'EOF'
unset PYTHONNOUSERSITE
EOF

echo "[INFO] Configured env to ignore ~/.local site-packages."

# ---------------- VERIFICATION ---------------- #
echo "[INFO] Verifying package versions..."
python -c "import sys, numpy, pandas; \
print('Python exe:', sys.executable); \
print('NumPy     :', numpy.__version__, '->', numpy.__file__); \
print('pandas    :', pandas.__version__, '->', pandas.__file__); \
assert numpy.__version__.startswith('1.22'), 'Unexpected NumPy version!'; \
assert pandas.__version__.startswith('2.0.2'), 'Unexpected pandas version!'"

echo "[INFO] Environment '$ENV_NAME' is ready."
echo "[INFO] Activate with: conda activate $ENV_NAME"
