GRP=NeglectedDiagnostics
PLAN=${GRP}Plan
ACR=neglecteddiagnosticsregistry
APP=neglected-diagnostics
LOC=westus3
IMG=$ACR.azurecr.io/neglected-diagnostics
# Should be run from the root of the project

if [[ "$SUBSCRIPTION" == "" ]]; then
       echo "Please set the SUBSCRIPTION environment variable to the name of your subscription"
       exit 1
fi

set -e
# Build a local wheel to help with installing dependencies in the container. See Dockerfile
rm -rf dist/*.whl
python -m build --wheel .
# Can be used to build/run locally:
# docker build -t ndiag --platform=linux/arm64 -f deployment/Dockerfile .
# docker run --rm -p 127.0.0.1:8501:80 ndiag

# ACR build is preferrably to docker build since it builds in azure and avoids having to build cross-platform locally
az account set --subscription "${SUBSCRIPTION}" # Use SSEC subscription
az group create -l westus3 -n $GRP
az acr create --name $ACR --resource-group $GRP --sku basic --admin-enabled true
az acr login --name $ACR
az acr build --registry $ACR --resource-group $GRP --image $APP -f deployment/Dockerfile .
az group create -l $LOC -n $GRP
az acr create --name $ACR --resource-group $GRP --sku basic --admin-enabled true
az appservice plan create -g $GRP -n $PLAN -l $LOC --is-linux --sku B2
az webapp create -g $GRP -p $PLAN -n $APP -i $IMG
az webapp deployment container config  -g $GRP -n $APP --enable-cd true