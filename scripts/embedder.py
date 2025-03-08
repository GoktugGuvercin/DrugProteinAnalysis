import re
import torch
from transformers import (
    T5Tokenizer,
    T5EncoderModel,
    BertTokenizer,
    BertModel,
)


class ProtT5Embedder:
    """It builds Prot T5 XL Uniref 50 Model for protein embeddings."""
    def __init__(self, device: str) -> None:
        self.device = device
        ckpt_name = "Rostlab/prot_t5_xl_half_uniref50-enc"
        self.tokenizer = T5Tokenizer.from_pretrained(ckpt_name, do_lower_case=False)
        self.model = T5EncoderModel.from_pretrained(ckpt_name).to(device)

        if self.device == "cpu":
            self.model.to(torch.float32)

    def compute_embeds(self, prot_seqs: list[str]) -> torch.Tensor:
        """It compute embeddings of protein sequences.
        Args:
          - prot_seqs: the amino acid sequences of proteins"""

        seqs = [" ".join(list(re.sub(r"[UZOB]", "X", seq))) for seq in prot_seqs]
        ids = self.tokenizer(seqs, padding="longest", return_tensors="pt")

        for k, v in ids.items():
            ids[k] = v.to(self.device)

        with torch.no_grad():
            output = self.model(**ids)

        embeds = []
        for i , seq in enumerate(prot_seqs):
            embed = output.last_hidden_state[i, :len(seq)].mean(dim=0)
            embeds.append(embed)

        return torch.stack(embeds)


class ProtTransEmbedder:
    """It builds ProtBert Model for protein embeddings."""
    def __init__(self, device: str) -> None:
        self.device = device
        ckpt_name = "Rostlab/prot_bert"
        self.tokenizer = BertTokenizer.from_pretrained(ckpt_name, do_lower_case=False)
        self.model = BertModel.from_pretrained(ckpt_name).to(device)

        if self.device == "cpu":
            self.model.to(torch.float32)

    def compute_res_embeds(self, prot_seqs: list[str]) -> torch.Tensor:
        """It compute mean residue embeddings of protein sequences.
        Args:
          - prot_seqs: the amino acid sequences of proteins"""

        seqs = [" ".join(list(re.sub(r"[UZOB]", "X", seq))) for seq in prot_seqs]
        ids = self.tokenizer(seqs, padding="longest", return_tensors="pt")

        for k, v in ids.items():
            ids[k] = v.to(self.device)

        with torch.no_grad():
            output = self.model(**ids)

        embeds = []
        for i , seq in enumerate(prot_seqs):
            embed = output.last_hidden_state[i, 1:len(seq) + 1].mean(dim=0)
            embeds.append(embed)

        return torch.stack(embeds)

    def get_cls_embeds(self, prot_seqs: list[str]) -> torch.Tensor:
        """It compute cls embeddings of protein sequences.
        Args:
          - prot_seqs: the amino acid sequences of proteins"""

        seqs = [" ".join(list(re.sub(r"[UZOB]", "X", seq))) for seq in prot_seqs]
        ids = self.tokenizer(seqs, padding="longest", return_tensors="pt")

        for k, v in ids.items():
            ids[k] = v.to(self.device)

        with torch.no_grad():
            output = self.model(**ids)

        return output.last_hidden_state[:, 0, :]
