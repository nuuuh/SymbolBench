# import os
# import openai
# from openai import AzureOpenAI 
# import base64

# endpoint = "https://egm-openai-1.openai.azure.com/" 
# model_name = "gpt-4o-mini" 
# deployment = "gpt-4o-mini"
# subscription_key = "c481f70b3fd14f78b293c554932bba3e"
# api_version = "2024-12-01-preview"


# class OpenAIModel(object):
#     def __init__(self, model_name, device, dtype, cache_dir=None, **kwargs):
#         self.model_name = model_name
#         self.device = device
#         self.dtype = dtype
#         self.api_key = subscription_key
#         assert self.api_key is not None, "API key not found."

#         self.client = AzureOpenAI( api_version=api_version, azure_endpoint=endpoint, api_key=subscription_key)

#         self.top_p = 0.9 if "top_p" not in kwargs else kwargs["top_p"]
#         self.temperature = 1.0 if "temperature" not in kwargs else kwargs["temperature"]
#         self.max_length = 1024 if "max_length" not in kwargs else kwargs["max_length"]
#         self.num_return_sequences = 1 if "num_return_sequences" not in kwargs else kwargs["num_return_sequences"]
#         self.seed = None if "seed" not in kwargs else kwargs["seed"]

#     def get_api_key(self):
#         if "OPENAI_API_KEY" in os.environ:
#             return os.environ["OPENAI_API_KEY"]
#         elif self.api_key_path is not None:
#             with open(self.api_key_path, "r") as f:
#                 return f.read().strip()
#         return None


#     def get_org_id(self):
#         if "OPENAI_ORG_ID" in os.environ:
#             return os.environ["OPENAI_ORG_ID"]
#         elif self.organization_id_path is not None:
#             with open(self.organization_id_path, "r") as f:
#                 return f.read().strip()

#         return None
    
#     def encode_image(self, image_path):
#         with open(image_path, "rb") as image_file:
#             return base64.b64encode(image_file.read()).decode('utf-8')
    
#     def get_messages(self, prompt, splits=["system", "user"], image_files=None):
#         messages = []
#         for split in splits:
#             start_tag = f"<{split}>"
#             end_tag = f"</{split}>"
#             if start_tag not in prompt or end_tag not in prompt:
#                 continue

#             start_idx = prompt.find(start_tag)
#             end_idx = prompt.find(end_tag)

#             messages.append({
#                 "role": split,
#                 "content": prompt[start_idx + len(start_tag):end_idx].strip()
#             })

#         if len(messages) == 0:
#             messages.append({
#                 "role": "user",
#                 "content": prompt
#             })

#         if image_files is not None:
#             user_index = next((i for i, item in enumerate(messages) if item["role"] == "user"), None)
#             user_msg = messages[user_index]
#             user_msg["content"] = [{"type": "text", "text": user_msg["content"]}]
#             for image_file in image_files:
#                 base64_image = self.encode_image(image_file)
#                 user_msg["content"].append({
#                     "type": "image_url",  # Updated type to 'image_url'
#                     "image_url": {'url':f"data:image/jpeg;base64,{base64_image}"}
#                 })
#         return messages
    
    
#     def generate(self, prompt, additional_prompt=None, image_files=None, temperature=None, max_new_tokens=None):
#         if temperature is None:
#             temperature = self.temperature
#         if max_new_tokens is None:
#             max_new_tokens = self.max_length

#         # import ipdb; ipdb.set_trace()
#         messages = self.get_messages(prompt, image_files=[additional_prompt] if additional_prompt is not None else None)
#         # import ipdb; ipdb.set_trace()
#         response = self.client.chat.completions.create(
#             messages=messages,
#             model=self.model_name,
#             temperature=temperature,
#             max_tokens=max_new_tokens,
#             top_p=self.top_p,
#             n=self.num_return_sequences,
#             seed=self.seed
#         )

#         if self.num_return_sequences==1:
#             return response.choices[0].message.content.strip()
#         else:
#             return [choice.message.content.strip() for choice in response.choices]
        

import os
import openai
import base64

class OpenAIModel(object):
    def __init__(self, model_name, device, dtype, cache_dir=None, **kwargs):
        self.model_name = model_name
        self.device = device
        self.dtype = dtype

        self.api_key_path = None if "api_key_path" not in kwargs else kwargs["api_key_path"]
        self.api_key = os.environ.get("OPENAI_API_KEY") # self.get_api_key()
        self.organization_id_path = None if "organization_id_path" not in kwargs else kwargs["organization_id_path"]
        self.organization_id = 'org-wyDIkK2j9B7M55gEHnkjc1Vd' # self.get_org_id()
        assert self.api_key is not None, "API key not found."

        self.client = openai.Client(
            api_key=self.api_key,
            organization=self.organization_id
        )

        self.top_p = 0.9 if "top_p" not in kwargs else kwargs["top_p"]
        self.temperature = 1.0 if "temperature" not in kwargs else kwargs["temperature"]
        self.max_length = 1024 if "max_length" not in kwargs else kwargs["max_length"]
        self.num_return_sequences = 1 if "num_return_sequences" not in kwargs else kwargs["num_return_sequences"]
        self.seed = None if "seed" not in kwargs else kwargs["seed"]

    def get_api_key(self):
        if "OPENAI_API_KEY" in os.environ:
            return os.environ["OPENAI_API_KEY"]
        elif self.api_key_path is not None:
            with open(self.api_key_path, "r") as f:
                return f.read().strip()
        return None


    def get_org_id(self):
        if "OPENAI_ORG_ID" in os.environ:
            return os.environ["OPENAI_ORG_ID"]
        elif self.organization_id_path is not None:
            with open(self.organization_id_path, "r") as f:
                return f.read().strip()

        return None
    
    def encode_image(self, image_path):
        with open(image_path, "rb") as image_file:
            return base64.b64encode(image_file.read()).decode('utf-8')
    
    def get_messages(self, prompt, splits=["system", "user"], image_files=None):
        messages = []
        for split in splits:
            start_tag = f"<{split}>"
            end_tag = f"</{split}>"
            if start_tag not in prompt or end_tag not in prompt:
                continue

            start_idx = prompt.find(start_tag)
            end_idx = prompt.find(end_tag)

            messages.append({
                "role": split,
                "content": prompt[start_idx + len(start_tag):end_idx].strip()
            })

        if len(messages) == 0:
            messages.append({
                "role": "user",
                "content": prompt
            })

        if image_files is not None:
            user_index = next((i for i, item in enumerate(messages) if item["role"] == "user"), None)
            user_msg = messages[user_index]
            user_msg["content"] = [{"type": "text", "text": user_msg["content"]}]
            for image_file in image_files:
                base64_image = self.encode_image(image_file)
                user_msg["content"].append({
                    "type": "image_url",
                    "image_url": {'url': f"data:image/jpeg;base64,{base64_image}"}
                    })
        return messages
    
    
    def generate(self, prompt, additional_prompt=None, image_files=None, temperature=None, max_new_tokens=None):
        if temperature is None:
            temperature = self.temperature
        if max_new_tokens is None:
            max_new_tokens = self.max_length
        # import ipdb; ipdb.set_trace()
        if additional_prompt is not None:
            additional_prompt = [additional_prompt] if isinstance(additional_prompt, str) else additional_prompt
        messages = self.get_messages(prompt, image_files=additional_prompt if additional_prompt is not None else None)
        
        try:
            response = self.client.chat.completions.create(
                messages=messages,
                model=self.model_name,
                temperature=temperature,
                max_tokens=max_new_tokens,
                top_p=self.top_p,
                n=self.num_return_sequences,
                seed=self.seed
            )

            if self.num_return_sequences==1:
                return response.choices[0].message.content.strip()
            else:
                return [choice.message.content.strip() for choice in response.choices]
        except Exception as e:
            # If OpenAI client fails (for example the model is not supported), try Together.ai as a fallback.
            try:
                from together import Together
            except Exception:
                # together not available; re-raise the original exception for visibility
                raise e

            try:
                tclient = Together()  # auth defaults to TOGETHER_API_KEY env var

                # Convert messages to plain-text form expected by Together
                conv_messages = []
                for m in messages:
                    content = m.get("content")
                    content_text = ""
                    if isinstance(content, list):
                        parts = []
                        for item in content:
                            if isinstance(item, str):
                                parts.append(item)
                            elif isinstance(item, dict) and item.get("type") == "text" and "text" in item:
                                parts.append(item["text"])
                            elif isinstance(item, dict) and "text" in item:
                                parts.append(item["text"])
                            # skip images or non-text content for the fallback
                        content_text = "\n".join(parts).strip()
                    elif isinstance(content, dict) and "text" in content:
                        content_text = content["text"].strip()
                    else:
                        content_text = str(content).strip()

                    conv_messages.append({"role": m.get("role", "user"), "content": content_text})

                response = tclient.chat.completions.create(
                    model=self.model_name,
                    messages=conv_messages,
                    temperature=temperature,
                    max_tokens=max_new_tokens,
                )

                # import ipdb; ipdb.set_trace()

                if self.num_return_sequences == 1:
                    return response.choices[0].message.content.strip()
                else:
                    return [choice.message.content.strip() for choice in response.choices]
            except Exception:
                # If Together also fails, re-raise the original OpenAI exception for debugging
                raise e
        