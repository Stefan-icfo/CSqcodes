import requests
import database
from qcodes import load_by_guid
from qcodes.dataset.plotting import plot_dataset
import base64

webhook_url = "https://icfo.webhook.office.com/webhookb2/4b7a2db6-80fb-4665-aa98-e5f3c767dfd0@f78a768a-22ae-4432-9eb4-55ce4b73c8c3/IncomingWebhook/88afbcf8f36348fcb8c827cea3da222a/a0d3c321-43f1-412b-894b-d2515a6e3261"


def convert_dataset_to_image_data(dataset):
    axis, _ = plot_dataset(dataset)
    axis[0].figure.savefig('tmp.jpeg')

    with open('tmp.jpeg', 'rb') as image:
        encoded = base64.b64encode(image.read())

    dataurl = f"data:image/jpeg;base64,{encoded.decode('utf-8')}"

    return dataurl


def adaptive_card(datasaver, log_params):
    payload = {
        "type": "message",
        "attachments": [
            {
                "contentType": "application/vnd.microsoft.card.adaptive",
                "contentUrl": None,
                "content": {
                    "$schema": "http://adaptivecards.io/schemas/adaptive-card.json",
                    "type": "AdaptiveCard",
                    "version": "1.2",
                    "body": [
                        {
                            "type": "TextBlock",
                            "text": f"The measurement with the id {datasaver.run_id}",
                            "size": "Large",
                            "wrap": True
                        },
                        {
                            "type": "Image",
                            "url": convert_dataset_to_image_data(datasaver.dataset),
                            "width": '1500px',
                        },

                    ]
                }
            }
        ]
    }

    for param in log_params:
        payload['attachments'][0]['content']['body'].append({
            "type": "TextBlock",
            "text": f"{param.label} -> {param()} {param.unit}",
            "wrap": True
        })

    response = requests.post(webhook_url, json=payload)
    print(response.content)
