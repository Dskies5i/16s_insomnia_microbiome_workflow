#!/usr/bin/env bash
# ==============================================================================
#   NAME: setup_and_run.sh
# SCRIPT MAESTRO DE DESPLIEGUE Y EJECUCIÓN (Linux Mint / macOS - Agnóstico)
# Proyecto: 16s_insomnia
# PERMISO (SIEMPRE LA PRIMERA VEZ): chmod +x 00_setup/setup_and_run.sh
# USAGE: bash 00_setup/setup_and_run.sh
# WORK DIRECTORY: cd progra26/16s_insomnia
# ==============================================================================

GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}==================================================${NC}"
echo -e "${BLUE}  🚀 INICIANDO CONFIGURACIÓN (16S INSOMNIA)${NC}"
echo -e "${BLUE}==================================================${NC}"

#-------------------------------------------------------------------------------
# 0. Dependencias del sistema operativo (Solo Linux)
OS_TYPE="$(uname)"
if [ "${OS_TYPE}" = "Linux" ]; then
    if ! command -v cargo &> /dev/null || ! pkg-config --exists zlib &> /dev/null; then
        echo -e "${YELLOW}⚠️ Librerías de desarrollo del sistema o Rust faltantes. Instalando vía apt...${NC}"
        sudo apt update && sudo apt install -y build-essential zlib1g-dev libxml2-dev libfontconfig1-dev libfreetype6-dev cargo rustc
    fi
elif [ "${OS_TYPE}" = "Darwin" ]; then
    echo -e "${GREEN}🍏 Sistema detectado: macOS. Usando entorno nativo Darwin.${NC}"
fi

#-------------------------------------------------------------------------------
# 1. Verificar si Conda/Mamba está instalado
if ! command -v conda &> /dev/null; then
    echo -e "${RED}❌ Error: Conda/Miniforge no está instalado o no está en el PATH.${NC}"
    exit 1
fi

PKG_MGR="conda"
if command -v mamba &> /dev/null; then
    PKG_MGR="mamba"
    echo -e "${GREEN}⚡ Usando Mamba para gestión rápida de entornos.${NC}"
fi

eval "$(conda shell.bash hook)"

#-------------------------------------------------------------------------------
# 2. Gestionar el entorno qiime2-amplicon-2024.5 (Python 3.9 aislado)
ENV_NAME="qiime2-amplicon-2024.5"
CONDA_BASE=$(conda info --base)
ENV_PATH="${CONDA_BASE}/envs/${ENV_NAME}"

check_python_version() {
    local target_env_path="$1"
    if [ -f "${target_env_path}/bin/python" ]; then
        local py_ver
        py_ver=$("${target_env_path}/bin/python" --version 2>&1 | awk '{print $2}')
        if [[ "$py_ver" =~ ^3\.9 ]]; then
            return 0
        fi
    fi
    return 1
}

if [ -d "$ENV_PATH" ] && check_python_version "$ENV_PATH"; then
    echo -e "${GREEN}✅ Entorno '${ENV_NAME}' detectado con Python 3.9. Activando...${NC}"
    conda activate "${ENV_NAME}"
else
    if [ -d "$ENV_PATH" ]; then
        echo -e "${YELLOW}⚠️ El entorno '${ENV_NAME}' existe pero tiene una versión incompatible. Recreándolo...${NC}"
        conda env remove -n "${ENV_NAME}" -y
    else
        echo -e "${YELLOW}⚠️ El entorno '${ENV_NAME}' no está instalado. Creándolo directamente...${NC}"
    fi
    
    # Creación directa mediante canales oficiales de QIIME 2 (Python 3.9 exigido)
    mamba create -n "${ENV_NAME}" python=3.9 -y
    conda activate "${ENV_NAME}"
    mamba install -c https://packages.qiime2.org/qiime2/2024.5/amplicon/released -c conda-forge -c bioconda -c defaults qiime2-amplicon=2024.5 -y
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}❌ Error crítico: Falló la creación del entorno con QIIME2.${NC}"
        exit 1
    fi
fi

# Asegurar que estamos usando el Rscript interno del entorno Conda activo
CONDA_R_SCRIPT="${CONDA_PREFIX}/bin/Rscript"
if [ ! -f "${CONDA_R_SCRIPT}" ]; then
    CONDA_R_SCRIPT="Rscript"
fi

R_LIB="${CONDA_PREFIX}/lib/R/library"

#-------------------------------------------------------------------------------
# 3. Verificar directorios clave del proyecto
echo -e "${BLUE}📁 Verificando estructura local de directorios...${NC}"
mkdir -p outputs/tables outputs/figures outputs/models outputs/logs \
         outputs/function/figures outputs/function/tables \
         outputs/boruta outputs/lefse outputs/maaslin2 \
         data/pre_processed/unified resources/classifiers

#-------------------------------------------------------------------------------
# 4. Asegurar paquetes críticos de R vía Binarios de Conda (Incluye todas las dependencias de microbiomeMarker)
echo -e "${BLUE}🔄 Instalando y asegurando dependencias en R (vía binarios seguros de Conda)...${NC}"

"${PKG_MGR}" install -c conda-forge -c bioconda \
    r-ggpubr r-yaml r-fs r-readr r-dplyr r-ggplot2 r-vegan \
    r-httr r-rsqlite r-desctools r-boruta r-cvxr \
    r-httpuv r-xml r-png r-haven r-systemfonts r-gdtools r-shiny \
    r-gridsvg r-latticeextra r-patchwork r-pheatmap r-rcolorbrewer r-circlize \
    bioconductor-phyloseq bioconductor-complexheatmap bioconductor-aldex2 \
    bioconductor-ggtree bioconductor-ancombc bioconductor-maaslin2 \
    bioconductor-microbiomemarker -y

# [FIX] Verificación defensiva física en disco para TODAS las dependencias críticas y secundarias
echo -e "${BLUE}🔍 Verificando existencia física de paquetes críticos en disco...${NC}"
for pkg in httr RSQLite DescTools CVXR ggtree microbiomemarker; do
    pkg_lower=$(echo "${pkg}" | tr '[:upper:]' '[:lower:]')
    if [ ! -d "${R_LIB}/${pkg}" ] && [ ! -d "${R_LIB}/${pkg_lower}" ]; then
        echo -e "${YELLOW}⚠️ El paquete ${pkg} no está físicamente en disco. Forzando reinstalación limpia...${NC}"
        if [[ "${pkg}" == "ggtree" || "${pkg}" == "microbiomemarker" ]]; then
            "${PKG_MGR}" install -c bioconda "bioconductor-${pkg_lower}" --force-reinstall -y
        else
            "${PKG_MGR}" install -c conda-forge "r-${pkg_lower}" --force-reinstall -y
        fi
    else
        echo -e "${GREEN}✅ ${pkg} verificado correctamente en disco.${NC}"
        if [ ! -d "${R_LIB}/${pkg}" ] && [ -d "${R_LIB}/${pkg_lower}" ]; then
            ln -s "${R_LIB}/${pkg_lower}" "${R_LIB}/${pkg}" 2>/dev/null || true
        fi
    fi
done

# Validación final limpia en R
"${CONDA_R_SCRIPT}" -e '
Sys.setenv(RENV_ACTIVATE_PROJECT = "FALSE")

core_pkgs <- c("phyloseq", "vegan", "ggplot2", "dplyr", "ggpubr", "Boruta",
               "readr", "yaml", "fs", "microbiomeMarker", "Maaslin2", "circlize",
               "httr", "RSQLite", "DescTools", "CVXR", "ggtree")
missing <- core_pkgs[!sapply(core_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  message("⚠️ Atención, faltan todavía estos paquetes: ", paste(missing, collapse = ", "))
} else {
  message("🎉 ¡Todos los paquetes esenciales y dependencias de microbiomeMarker están listos y cargan sin errores!")
}
'

echo -e "${GREEN}✅ Configuración completada con éxito.${NC}"

#-------------------------------------------------------------------------------
# 5. Menú de ejecución
echo -e "\n${BLUE}==================================================${NC}"
echo -e "${BLUE}  ⚙️  ¿Qué deseas ejecutar?${NC}"
echo -e "  1) Ver guía / Recordatorio de Preprocesamiento (Manual)"
echo -e "  2) Ejecutar Metaanálisis y Análisis Funcional completo (R)"
echo -e "${BLUE}==================================================${NC}"
read -p "Selecciona una opción [1-2]: " option

case $option in
    1)
        echo -e "${BLUE}ℹ️ Módulo de Preprocesamiento (0_preprocessing):${NC}"
        echo -e "Recuerda que el preprocesamiento se ejecuta de forma manual cohorte por cohorte."
        ;;
    2)
        echo -e "${YELLOW}💡 Recordatorio al agregar nuevos datasets:${NC}"
        echo -e "Si has añadido nuevos datasets, recuerda limpiar los outputs previos"
        echo -e "para evitar que los checkpoints usen datos obsoletos."
        echo -e "  Ejecuta manualmente si lo necesitas:"
        echo -e "  ${RED}rm -f outputs/tables/merged_KO_metagenome.tsv${NC}"
        echo -e "  ${RED}rm -rf picrust2_out${NC}"
        echo -e "  ${RED}rm -rf outputs/tables/* outputs/figures/* outputs/models/* outputs/logs/* outputs/boruta/* outputs/lefse/* outputs/maaslin2/* outputs/function/* outputs/function_contrib/* outputs/dbrda/*${NC}\n"
        
        read -p "¿Deseas continuar con la ejecución del pipeline? (s/n): " confirm
        if [[ "${confirm}" =~ [sS] ]]; then
            echo -e "${BLUE}▶️ Ejecutando pipeline completo de metaanálisis (R)...${NC}"
            export RENV_ACTIVATE_PROJECT="FALSE"
            "${CONDA_R_SCRIPT}" 1_meta_analysis/run_all.R
            echo -e "${GREEN}🎉 ¡Metaanálisis y análisis funcional completados con éxito!${NC}"
        else
            echo -e "${YELLOW}Ejecución cancelada por el usuario.${NC}"
        fi
        ;;
    *)
        echo -e "${RED}Opción inválida. Saliendo.${NC}"
        exit 1
        ;;
esac