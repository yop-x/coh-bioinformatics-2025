import streamlit as st # type: ignore
import networkx as nx
import matplotlib.pyplot as plt
import tempfile
from coexp_update import get_coexp_network_updated  # type: ignore
import requests 


st.set_page_config(page_title="Co-Expression Explorer", layout="wide")
st.title("🧬 RNA–DNA Co-Expression Explorer")


st.markdown("""
This tool visualizes co-expression networks using iMARGI datasets.
Enter a DNA coordinate or gene name to see linked RNAs or DNAs.
""")


imargi_file = st.file_uploader("upload your 10/15-column iMARGI .bedpe.gz file", type=["gz", "bedpe"])


freq = st.slider("Minimum interaction count (frequency)", min_value=3, max_value=50, value=5)
query_input = st.text_input("Enter DNA gene symbol (e.g. LINC00607) or coordinate (e.g. chr1:1000000-1010000):")


if imargi_file is not None and query_input:
    st.success("Uploading file to Colab server...")

    # Send file and inputs to Colab FastAPI
    files = {"file": imargi_file}
    data = {
        "query": query_input,
        "freq": freq
    }

    # 🔁 Paste your public URL here!
    colab_url = "http://763c865d10af.ngrok-free.app/run"

    try:
        response = requests.post(colab_url, files=files, data=data)
        if response.status_code == 200:
            st.image(response.content, caption="Co-expression Network")
        else:
            st.error(f"Server error {response.status_code}: {response.text}")
    except Exception as e:
        st.error(f"Error sending request: {e}")