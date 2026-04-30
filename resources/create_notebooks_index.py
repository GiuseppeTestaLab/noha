import os
import argparse
import yaml
import markdown
from collections import defaultdict

INDEX_TEMPLATE_HEADER = """
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="X-UA-Compatible" content="ie=edge">
    <title>{short_title}</title>
    <link rel="stylesheet" href="./style.css">
    <link rel="icon" href="./favicon.ico" type="image/x-icon">
  </head>
  <body>
    <main>
      <header class="site-header">
        <h1>{title}</h1>
        <div class="subtitle">{subtitle}</div>
      </header>
"""

SECTION_TEMPLATE = """
      <section class="notebook-section">
        <h2 class="section-title">{section}</h2>
        <div class="notebook-grid">
{items}
        </div>
      </section>
"""

CARD_TEMPLATE = """
          <article class="notebook-card">
            <h3><a href="{href}">{title}</a></h3>
            <div class="description">{description}</div>
          </article>
"""

INDEX_TEMPLATE_FINAL = """
      <footer class="site-footer">
        {final}
      </footer>
    </main>
    <script src="index.js"></script>
  </body>
</html>
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("configyaml")
    args = parser.parse_args()

    with open(args.configyaml) as f:
        config_dict = yaml.load(f, Loader=yaml.SafeLoader)

    html_subtitle = markdown.markdown(config_dict.get("subtitle", ""))
    html_final = markdown.markdown(config_dict.get("final", ""))

    print(INDEX_TEMPLATE_HEADER.format(
        title=config_dict["title"],
        short_title=config_dict["short_title"],
        subtitle=html_subtitle
    ))

    sections = defaultdict(list)

    for key, page in config_dict["pages"].items():
        section = page.get("section", "Other")
        sections[section].append((key, page))

    for section, items in sections.items():
        cards = []

        for key, page in items:
            html_description = markdown.markdown(page["description"])
            path_to_file = page["notebook"].replace(".ipynb", ".html")
            file_name = os.path.basename(path_to_file)

            cards.append(CARD_TEMPLATE.format(
                title=key,
                href=file_name,
                description=html_description
            ))

        print(SECTION_TEMPLATE.format(
            section=section,
            items="".join(cards)
        ))

    print(INDEX_TEMPLATE_FINAL.format(final=html_final))

if __name__ == '__main__':
    main()
