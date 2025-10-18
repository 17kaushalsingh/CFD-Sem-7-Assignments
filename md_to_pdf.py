import markdown

# Read markdown file
with open('report.md', 'r', encoding='utf-8') as f:
    md_content = f.read()

# Convert to HTML
html = markdown.markdown(md_content, extensions=['tables', 'fenced_code'])

# Create full HTML document with styling
html_doc = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            line-height: 1.6;
            max-width: 900px;
            margin: 40px auto;
            padding: 0 20px;
            color: #333;
        }}
        h1 {{ border-bottom: 2px solid #333; padding-bottom: 10px; }}
        h2 {{ border-bottom: 1px solid #ccc; padding-bottom: 8px; margin-top: 30px; }}
        h3 {{ margin-top: 25px; }}
        code {{
            background: #f4f4f4;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Courier New', monospace;
        }}
        pre {{
            background: #f4f4f4;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }}
        th {{
            background-color: #f4f4f4;
            font-weight: bold;
        }}
        hr {{ border: 0; border-top: 1px solid #ddd; margin: 30px 0; }}
    </style>
</head>
<body>
{html}
</body>
</html>
"""

# Write to HTML file
with open('documentation.html', 'w', encoding='utf-8') as f:
    f.write(html_doc)

print("✓ Created documentation.html")
print("Open it in your browser and use Print -> Save as PDF")