# gex_predictor/__init__.py
import ssl
import certifi
import os

os.environ["REQUESTS_CA_BUNDLE"] = certifi.where()
os.environ["SSL_CERT_FILE"] = certifi.where()
ssl._create_default_https_context = ssl.create_default_context(
    cafile=certifi.where()
)

from .GEx_Predictor import GEx_Predictor
